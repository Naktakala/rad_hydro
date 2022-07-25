#include "rhsolverA.h"

#include "ChiMath/chi_math_banded_solvers.h"

#include "ChiTimer/chi_timer.h"

void chi_radhydro::SolverA_GDCN::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  //======================================== Get options
  const auto num_timesteps     = basic_options("max_timesteps").IntegerValue();
  const auto delta_t_max       = basic_options("maximum_dt"   ).FloatValue();
  const auto CFL               = basic_options("CFL"          ).FloatValue();
  const auto max_time          = basic_options("max_time"     ).FloatValue();
  const auto output_times_str  = basic_options("export_times" ).StringValue();

  std::vector<double> output_times = MakeOutputTimesFromStr(output_times_str);

  const auto& CCL   = m_scalar_fields.at("cell_char_length");

  m_kappa_s_function = basic_options("kappa_s_function").StringValue();
  m_kappa_a_function = basic_options("kappa_a_function").StringValue();

  //======================================== Allocate work data
  const size_t num_local_nodes = m_grid->local_cells.size();
  std::vector<UVector>  U_n(num_local_nodes);
  std::vector<UVector>  U_n_star(num_local_nodes);
  std::vector<UVector>  U_nph(num_local_nodes);
  std::vector<UVector>  U_np1(num_local_nodes);

  VecDbl                rad_E_n(num_local_nodes);
  VecDbl                rad_E_n_star(num_local_nodes);
  VecDbl                rad_E_nph(num_local_nodes);
  VecDbl                rad_E_np1(num_local_nodes);

  VecDbl                kappa_a_n(num_local_nodes, 0.0);
  VecDbl                kappa_t_n(num_local_nodes, 0.0);

  VecDbl                kappa_a_nph(num_local_nodes, 0.0);
  VecDbl                kappa_t_nph(num_local_nodes, 0.0);

  //======================================== Develop U and radE from field funcs
  FieldsToU(m_scalar_fields.at("rho"),
            m_scalar_fields.at("u"),
            m_scalar_fields.at("v"),
            m_scalar_fields.at("w"),
            m_scalar_fields.at("e"),
            U_n);

  rad_E_n = m_scalar_fields.at("radE");

  //======================================== Diagnostics on initial conditions
  auto PrintVal = [](const char* formatting, double value)
  {
    char buff[100];
    sprintf(buff, formatting, value);
    return std::string(buff);
  };
  auto PrintUradE = [this,PrintVal](const UVector& U, const double radE)
  {
    const double rho = U[RHO];
    const double u   = chi_radhydro::VelocityFromCellU(U).Norm();
    const double e   = chi_radhydro::InternalEnergyFromCellU(U);
    const double T   = e/m_Cv;

    std::stringstream ostr;
    ostr << PrintVal(" rho  %8.8e", rho);
    ostr << PrintVal(" u    %8.8e", u);
    ostr << PrintVal(" T    %8.8e", T);
    ostr << PrintVal(" e    %8.8e", e);
    ostr << PrintVal(" radE %8.8e", radE);
    ostr << PrintVal(" E    %8.8e", U[MAT_E]);

    return ostr.str();
  };
  {
    auto p_L = chi_radhydro::IdealGasPressureFromCellU(U_n.front(), m_gamma);
    auto p_R = chi_radhydro::IdealGasPressureFromCellU(U_n.back(), m_gamma);

    auto u_L = chi_radhydro::VelocityFromCellU(U_n.front()).Norm();
    auto u_R = chi_radhydro::VelocityFromCellU(U_n.back()).Norm();

    auto radEL = m_scalar_fields["radE"].front();
    auto radER = m_scalar_fields["radE"].back();

    auto F_L = chi_radhydro::MakeFWithRadE(U_n.front(), p_L, radEL, Vec3(0,0,1));
    auto F_R = chi_radhydro::MakeFWithRadE(U_n.back() , p_R, radER, Vec3(0,0,1));

    auto F_net = F_R - F_L;
    double Frad_E_net = radEL*u_L - radER*u_R;

    chi::log.Log() << "radER_def: " << PrintVal("%8.8e",radEL*u_L/u_R);

    chi::log.Log() << "F_L : " << F_L.PrintS() << " " << PrintVal("%8.8e",radEL);
    chi::log.Log() << "F_R : " << F_R.PrintS() << " " << PrintVal("%8.8e",radER);
    chi::log.Log() << "Fnet: " << F_net.PrintS() << " " << Frad_E_net;
    const size_t N = num_local_nodes;
    std::stringstream ostr;
    ostr << PrintVal(" rho0  %8.8e", m_scalar_fields["rho"        ][0]);
    ostr << PrintVal(" u0    %8.8e", m_scalar_fields["w"          ][0]);
    ostr << PrintVal(" T0    %8.8e", m_scalar_fields["temperature"][0]);
    ostr << PrintVal(" e0    %8.8e", m_scalar_fields["e"          ][0]);
    ostr << PrintVal(" radE0 %8.8e", m_scalar_fields["radE"       ][0]);
    ostr << PrintVal(" E0    %8.8e", U_n[0][MAT_E]);
    ostr << "\n";
    ostr << PrintVal(" rho1  %8.8e", m_scalar_fields["rho"        ][N - 1]);
    ostr << PrintVal(" u1    %8.8e", m_scalar_fields["w"          ][N - 1]);
    ostr << PrintVal(" T1    %8.8e", m_scalar_fields["temperature"][N - 1]);
    ostr << PrintVal(" e1    %8.8e", m_scalar_fields["e"          ][N - 1]);
    ostr << PrintVal(" radE1 %8.8e", m_scalar_fields["radE"       ][N - 1]);
    ostr << PrintVal(" E1    %8.8e", U_n[N - 1][MAT_E]);

    chi::log.Log() << ostr.str();
  }

  //======================================== Init shock location tracker
  double shock_ref_speed = 0.0;
  for (double uval : m_scalar_fields.at("w"))
    shock_ref_speed = std::max(shock_ref_speed, uval);

  double shock_location = 0.5;
  typedef std::pair<double,double> TimePosPair;
  std::vector<TimePosPair> shock_locations_at_output_times;

  //======================================== Initialize system energy
  double system_energy_change = 0.0;
  SystemEnergy new_system_energy = ComputeSysEnergyChange(0.0, U_n, rad_E_n,
                                                          {}, {},
                                                          m_bc_settings);
  SystemEnergy old_system_energy;

  //======================================== Pickup references
  SimRefs sim_refs{*m_grid,
                   *m_fv,
                   m_bc_settings,
                   m_Cv, m_gamma,
                   m_kappa_s_function, m_kappa_a_function};

  //######################################## Start iterations
  chi_objects::ChiTimer timer;
  std::unordered_map<std::string, double> timing_info;


  double time = 0.0; //Total time
  int n=0;           //Iteration count
  while (true)
  {
    timer.Reset();
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, m_gamma, CCL, CFL);
    const double dt       = std::min(delta_t_max, dt_hydro);

    timing_info["delta_t"] += GetResetTimer(timer);

    //================================= Print iteration info
    {
      using namespace std;
      chi::log.Log() << "Timestep " << setw(7) << n
                     << " with dt=" << setw(7)
                     << scientific << setprecision(1)<< dt
                     << " t=" << setw(7) << fixed << setprecision(4) << time
                     << " shock_loc: " << shock_location
                     << " delta_energy: " << setw(12)
                     << scientific << setprecision(4)
                     << system_energy_change
        << " matE_adv: " << setw(12)
        << scientific << setprecision(8) << new_system_energy.me_adv
        << " radE_adv: " << setw(12)
        << scientific << setprecision(8) << new_system_energy.re_adv;
    }
    old_system_energy = new_system_energy;

    //================================= Compute gradients at n
    //                                  With slope limiters
    const auto grad_U_n     = ComputeUGradients(U_n, sim_refs);
    const auto grad_rad_E_n = ComputeRadEGradients(rad_E_n, sim_refs);
    timing_info["gradient_n"] += GetResetTimer(timer);

    //================================= Compute kappas at n
    chi_radhydro::ComputeCellKappas(*m_grid, m_Cv, U_n,
                                    m_kappa_s_function, m_kappa_a_function,
                                    kappa_a_n, kappa_t_n);

    timing_info["compute_kappa_n"] += GetResetTimer(timer);

    //================================= Predictor
    Predictor(sim_refs,
              kappa_a_n, kappa_t_n,
              dt,
              U_n, grad_U_n, U_nph,
              rad_E_n, grad_rad_E_n, rad_E_nph);

//    //=================================== Advection of U and rad_E
//    MHM_HydroRadEPredictor(    //Defined in chi_radhydro
//      sim_refs,                //Stuff
//      /*tau=*/1/(0.5*dt),      //tau input
//      U_n, grad_U_n,           //Hydro inputs
//      rad_E_n, grad_rad_E_n,   //RadE inputs
//      U_n_star, rad_E_n_star   //Outputs
//    );
//
//    //=================================== Update density and momentum to nph
//    DensityMomentumUpdateWithRadMom(            //Defined in chi_radhydro
//      sim_refs,                                 //Stuff
//      kappa_t_n,                                //Kappa_t for interpolation
//      /*tau=*/1/(0.5*dt),                       //tau input
//      U_n, U_n_star,                            //Hydro inputs
//      rad_E_n,                                  //RadE input
//      U_nph                                     //Output
//    );
//
//    //=================================== Internal energy and radiation energy
//    {
//      std::vector<double> k5,k6;
//      MatDbl A;
//      VecDbl b;
//
//      AssembleGeneralEnergySystem(
//        sim_refs,
//        kappa_a_n, kappa_t_n,
//        kappa_a_n, kappa_t_n,                               //Stuff
//        /*tau=*/1/(0.5*dt), /*theta1=*/1.0, /*theta2=*/0.0, //tau, theta input
//        U_n, U_n, U_n_star, U_nph, grad_U_n,                //Hydro inputs
//        rad_E_n, rad_E_n, rad_E_n_star, grad_rad_E_n,       //RadE inputs
//        k5, k6, A, b);                                      //Outputs
//
//      rad_E_nph = chi_math::TDMA(A,b);
//
//      //============================ Back sub for internal energy
//      for (const auto& cell_c : m_grid->local_cells)
//      {
//        const uint64_t c     = cell_c.local_id;
//        double rho_c_nph     = U_nph[c][RHO];
//        Vec3   u_c_nph       = VelocityFromCellU(U_nph[c]);
//        double u_abs_sqr_nph = u_c_nph.NormSquare();
//        double e_c_nph       = k5[c]*rad_E_nph[c] + k6[c];
//
//        double& E_c_nph = U_nph[c](MAT_E);
//
//        E_c_nph = rho_c_nph*(0.5 * u_abs_sqr_nph + e_c_nph);
//      }//for cell c
//    }//scope internal e and rad_E

    timing_info["predictor"] += GetResetTimer(timer);

    //================================= Compute gradients at nph
    // With slope limiters
    const auto& grad_U_nph     = grad_U_n;
    const auto& grad_rad_E_nph = grad_rad_E_n;

    timing_info["gradient_nph"] += GetResetTimer(timer);

    //=================================== Compute kappas at nph
    chi_radhydro::ComputeCellKappas(*m_grid, m_Cv, U_nph,
                                    m_kappa_s_function, m_kappa_a_function,
                                    kappa_a_nph, kappa_t_nph);

    timing_info["compute_kappa_nph"] += GetResetTimer(timer);

    //================================= Corrector
    CorrectorHydroAndMom(sim_refs,
                         kappa_a_n, kappa_t_n,
                         kappa_a_nph, kappa_t_nph,
                         dt,
                         U_n, U_nph, grad_U_nph, U_np1,
                         rad_E_n, rad_E_nph, grad_rad_E_nph, rad_E_np1);

    timing_info["corrector"] += GetResetTimer(timer);


    //================================= Compute energy conservation
    new_system_energy = ComputeSysEnergyChange(dt, U_n, rad_E_n,
                                               grad_U_n, grad_rad_E_n,
                                               m_bc_settings);
    system_energy_change =
      new_system_energy.Emat - old_system_energy.Emat +
      new_system_energy.Erad - old_system_energy.Erad +
      new_system_energy.me_adv +
      new_system_energy.re_adv;
    old_system_energy.Emat += new_system_energy.me_adv;
    old_system_energy.Erad += new_system_energy.re_adv;

    //================================= Copy new solution to old
    U_n = U_np1;
    rad_E_n = rad_E_np1;
    time += dt;

    //================================= Finding shock location
    {
      size_t k=0;
      for (const auto& Uval : U_n)
      {
        if ((Uval[3]/Uval[0]) < (0.95*shock_ref_speed))
        {
          shock_location = m_grid->local_cells[k].centroid.z;
          break;
        }
        ++k;
      }
    }//scope: Finding shock location

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

          shock_locations_at_output_times.emplace_back(time, shock_location);

          ProcessAndExportFields(file_name,U_n,rad_E_n,num_local_nodes,
                                 {{"time",time},
                                  {"e_balance",system_energy_change}});
        }
      }//for output_time
    }

    ++n;

    if (num_timesteps > 0)
      if (n >= num_timesteps) break;
    if (max_time > 0.0)
      if (time >= max_time) break;
  }//for n
  //######################################## Done iterating

  chi::log.Log() << "Outflow conditions: \n"
                 << PrintUradE(U_n.back(), rad_E_n.back());

  //======================================== Print timing info
  {
    std::stringstream ostr;
    ostr << "Timing info:\n";
    for (const auto& timing_str_info : timing_info)
    {
      char buff1[100], buff2[100];
      sprintf(buff1, "%20s", timing_str_info.first.c_str());
      sprintf(buff2, "%20.2f", timing_str_info.second);
      ostr << buff1 << " " << buff2;
      ostr << "\n";
    }
    chi::log.Log() << ostr.str();
  }

  //======================================== Print shock locations
  {
    using namespace std;
    stringstream ostr;
    ostr << "Shock locations:\n";
    for (const auto& time_pos_pair : shock_locations_at_output_times)
      ostr << "t=" << setw(7) << fixed << setprecision(4) << time_pos_pair.first
           << " x="
           << setw(11) << fixed << setprecision(8) << time_pos_pair.second
           << "\n";
    if (shock_locations_at_output_times.size() > 3)
    {
      const auto t_u_i = shock_locations_at_output_times[1];
      const auto t_u_f = shock_locations_at_output_times.back();
      const double dt = t_u_f.first - t_u_i.first;
      const double du = t_u_f.second - t_u_i.second;

      ostr << "Estimated shock speed: "
           << setw(11) << fixed << setprecision(8) << du/dt;
    }
    chi::log.Log() << ostr.str();
  }

  //======================================== Parse U and radE to field-functions
  UToFields(U_n,
            m_scalar_fields.at("rho"),
            m_scalar_fields.at("u"),
            m_scalar_fields.at("v"),
            m_scalar_fields.at("w"),
            m_scalar_fields.at("e"),
            m_scalar_fields.at("p"),
            m_gamma);

  m_scalar_fields.at("radE") = rad_E_n;

  //======================================== Compute material temperature
  auto& temperature = m_scalar_fields.at("temperature");

  const auto& e  = m_scalar_fields.at("e");

  for (uint64_t c=0; c < num_local_nodes; ++c)
    temperature[c] = e[c] / m_Cv;

  //======================================== Print outputs
  PrintRawOutput("ZRawOutput.txt",{});

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}