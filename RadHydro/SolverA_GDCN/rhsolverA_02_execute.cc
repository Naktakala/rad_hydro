#include "rhsolverA.h"

#include "ChiTimer/chi_timer.h"

void chi_radhydro::SolverA_GDCN::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  //======================================== Define lambdas
  auto ProcessAndExportFields = [this](const std::string& file_name,
                                       const std::vector<UVector>& U,
                                       const std::vector<double>& rad_E,
                                       const std::vector<double>& Cv,
                                       const size_t num_local_nodes)
  {
    //======================================== Parse U and radE to field-functions
    UToFields(U,
              scalar_fields.at("rho"),
              scalar_fields.at("u"),
              scalar_fields.at("v"),
              scalar_fields.at("w"),
              scalar_fields.at("e"),
              scalar_fields.at("p"),
              scalar_fields.at("gamma"));

    scalar_fields.at("radE") = rad_E;

    //======================================== Compute material temperature
    auto& temperature = scalar_fields.at("temperature");

    const auto& e  = scalar_fields.at("e");

    for (uint64_t c=0; c < num_local_nodes; ++c)
      temperature[c] = e[c] / Cv[c];

    //======================================== Print outputs
    PrintRawOutput(file_name);
  };

  //======================================== Allocate work data
  const size_t num_local_nodes = grid->local_cells.size();
  std::vector<UVector>  U_n(num_local_nodes);
  std::vector<UVector>  U_nph(num_local_nodes);
  std::vector<UVector>  U_np1(num_local_nodes);

  VecDbl                rad_E_n(num_local_nodes);
  VecDbl                rad_E_nph(num_local_nodes);
  VecDbl                rad_E_np1(num_local_nodes);

  VecDbl                kappa_a_n(num_local_nodes, 0.0);
  VecDbl                kappa_t_n(num_local_nodes, 0.0);

  VecDbl                kappa_a_nph(num_local_nodes, 0.0);
  VecDbl                kappa_t_nph(num_local_nodes, 0.0);

  //======================================== Develop U and radE from field funcs
  FieldsToU(scalar_fields.at("rho"),
            scalar_fields.at("u"),
            scalar_fields.at("v"),
            scalar_fields.at("w"),
            scalar_fields.at("e"),
            U_n);

  rad_E_n = scalar_fields.at("radE");

  {
    auto PrintF = [](const char* formatting, double value)
    {
      char buff[100];
      sprintf(buff, formatting, value);
      return std::string(buff);
    };
    double gamma = scalar_fields.at("gamma")[0];
    auto p_L = chi_radhydro::IdealGasPressureFromCellU(U_n.front(), gamma);
    auto p_R = chi_radhydro::IdealGasPressureFromCellU(U_n.back(), gamma);

    auto radEL = scalar_fields["radE"].front();
    auto radER = scalar_fields["radE"].back();

    auto F_L = chi_radhydro::MakeFWithRadE(U_n.front(), p_L, radEL, Vec3(0,0,1));
    auto F_R = chi_radhydro::MakeFWithRadE(U_n.back() , p_R, radER, Vec3(0,0,1));

    auto F_net = F_R - F_L;

    chi::log.Log() << "F_L : " << F_L.PrintS();
    chi::log.Log() << "F_R : " << F_R.PrintS();
    chi::log.Log() << "Fnet: " << F_net.PrintS();
    const size_t N = num_local_nodes;
    std::stringstream ostr;
    ostr << PrintF(" rho0  %8.8e", scalar_fields["rho"        ][0]);
    ostr << PrintF(" u0    %8.8e", scalar_fields["w"          ][0]);
    ostr << PrintF(" T0    %8.8e", scalar_fields["temperature"][0]);
    ostr << PrintF(" e0    %8.8e", scalar_fields["e"          ][0]);
    ostr << PrintF(" radE0 %8.8e", scalar_fields["radE"       ][0]);
    ostr << "\n";
    ostr << PrintF(" rho1  %8.8e", scalar_fields["rho"        ][N-1]);
    ostr << PrintF(" u1    %8.8e", scalar_fields["w"          ][N-1]);
    ostr << PrintF(" T1    %8.8e", scalar_fields["temperature"][N-1]);
    ostr << PrintF(" e1    %8.8e", scalar_fields["e"          ][N-1]);
    ostr << PrintF(" radE1 %8.8e", scalar_fields["radE"       ][N-1]);

    chi::log.Log() << ostr.str();

  }
//  return;

  //======================================== Get options
  const auto num_timesteps     = basic_options("max_timesteps").IntegerValue();
  const auto delta_t_max       = basic_options("maximum_dt"   ).FloatValue();
  const auto CFL               = basic_options("CFL"          ).FloatValue();
  const auto max_time          = basic_options("max_time"     ).FloatValue();
  const auto output_times_str  = basic_options("export_times" ).StringValue();

  std::vector<double> output_times;
  {
    std::istringstream iss(output_times_str);
    while (iss.tellg() >= 0)
    {
      double value; iss >> value;
      if (not (iss.rdstate() & std::istringstream::failbit))
        output_times.push_back(value);
    }
    std::stringstream oss;
    for (double val : output_times)
      oss << val << " ";
    chi::log.Log() << "Requested output times: " << oss.str();
  }


  const auto& CCL   = scalar_fields.at("cell_char_length");
  const auto& gamma = scalar_fields.at("gamma");
  const auto& Cv    = scalar_fields.at("Cv");

  auto kappa_s_function = basic_options("kappa_s_function").StringValue();
  auto kappa_a_function = basic_options("kappa_a_function").StringValue();

  //======================================== Init shock location tracker
  double shock_ref_speed = 0.0;
  for (double uval : scalar_fields.at("w"))
    shock_ref_speed = std::max(shock_ref_speed, uval);

  double shock_location = 0.5;

  //######################################## Start iterations
  auto GetResetTimer = [](chi_objects::ChiTimer& timer)
  {
    double time = timer.GetTime();
    timer.Reset();
    return time;
  };
  chi_objects::ChiTimer timer;
  std::map<std::string, double> timing_info;
  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    timer.Reset();
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, gamma, CCL, CFL);
    const double dt       = std::min(delta_t_max, dt_hydro);

    timing_info["delta_t"] += GetResetTimer(timer);

    chi::log.Log() << "Timestep " << n
                   << " with dt=" << dt
                   << " t=" << time
                   << " shock_loc: " << shock_location;

    //================================= Compute gradients at n
    // With slope limiters
    auto grad_U_n     = ComputeUGradients(U_n,*grid,*fv,bc_settings);
    auto grad_rad_E_n = ComputeRadEGradients(rad_E_n,*grid,*fv,bc_settings);
    timing_info["gradient_n"] += GetResetTimer(timer);

    //================================= Compute kappas at n
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_n,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_n,kappa_t_n);

    timing_info["compute_kappa_n"] += GetResetTimer(timer);

    //================================= Predictor
    Predictor(bc_settings, kappa_a_n, kappa_t_n,
              dt,
              U_n,grad_U_n,U_nph,
              rad_E_n,grad_rad_E_n,rad_E_nph);

    timing_info["predictor"] += GetResetTimer(timer);

    //================================= Compute gradients at nph
    // With slope limiters
    auto grad_U_nph     = ComputeUGradients(U_nph,*grid,*fv,bc_settings);
    auto grad_rad_E_nph = ComputeRadEGradients(rad_E_nph,*grid,*fv,bc_settings);

    timing_info["gradient_nph"] += GetResetTimer(timer);

    //=================================== Compute kappas at nph
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_nph,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_nph,kappa_t_nph);

    timing_info["compute_kappa_nph"] += GetResetTimer(timer);

    //================================= Corrector
    Corrector(bc_settings, kappa_a_n, kappa_t_n,
              kappa_a_nph, kappa_t_nph,
              dt,
              U_n,U_nph,grad_U_nph,U_np1,
              rad_E_n,rad_E_nph,grad_rad_E_nph,rad_E_np1);

    timing_info["corrector"] += GetResetTimer(timer);

    //================================= Copy new solution to old
    U_n = U_np1;
    rad_E_n = rad_E_np1;
    time += dt;

    //================================= Processing output times
    {
      int t = 0;
      for (double output_time : output_times)
      {
        ++t;
        char buffer[100];
        sprintf(buffer, "%03d",t);
        if ((time-dt) < output_time and time >= output_time)
        {
          const auto file_prefix = basic_options("output_prefix").StringValue();
          const auto file_name = file_prefix + "t" + buffer + ".txt";

          ProcessAndExportFields(file_name,U_n,rad_E_n,Cv,num_local_nodes);
        }
      }//for output_time
    }

    //================================= Finding shock location
    {
      size_t k=0;
      for (const auto& Uval : U_n)
      {
        ++k;
        if ((Uval[3]/Uval[0]) < (0.95*shock_ref_speed))
        {
          shock_location = double(k)/double(num_local_nodes);
          break;
        }
      }
    }//scope: Finding shock location

    if (max_time > 0.0)
      if (time >= max_time) break;
  }//for n
  //######################################## Done iterating

  //======================================== Print timing info
  {
    std::stringstream ostr;
    ostr << "Timing info:\n";
    for (const auto& timing_str_info : timing_info)
    {
      ostr << timing_str_info.first << " " << timing_str_info.second;
      ostr << "\n";
    }
    chi::log.Log() << ostr.str();
  }

  //======================================== Parse U and radE to field-functions
  UToFields(U_n,
            scalar_fields.at("rho"),
            scalar_fields.at("u"),
            scalar_fields.at("v"),
            scalar_fields.at("w"),
            scalar_fields.at("e"),
            scalar_fields.at("p"),
            scalar_fields.at("gamma"));

  scalar_fields.at("radE") = rad_E_n;

  //======================================== Compute material temperature
  auto& temperature = scalar_fields.at("temperature");

  const auto& e  = scalar_fields.at("e");

  for (uint64_t c=0; c < num_local_nodes; ++c)
    temperature[c] = e[c] / Cv[c];

  //======================================== Print outputs
  PrintRawOutput("ZRawOutput.txt");

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}