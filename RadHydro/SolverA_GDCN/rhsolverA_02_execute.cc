#include "rhsolverA.h"

void chi_radhydro::SolverA_GDCN::Execute()
{
  auto PrintField = [this](const std::string& field_name)
  {
    std::cout << field_name << ": ";
    const auto& field = scalar_fields.at(field_name);
    for (auto val : field) std::cout << val << " ";
    std::cout << "\n";
  };

  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  const size_t num_nodes_local = grid->local_cells.size();
  std::vector<UVector>  U_n(num_nodes_local);
  std::vector<UVector>  U_nph(num_nodes_local);
  std::vector<UVector>  U_np1(num_nodes_local);

  std::vector<double>   rad_E_n(num_nodes_local);
  std::vector<double>   rad_E_nph(num_nodes_local);
  std::vector<double>   rad_E_np1(num_nodes_local);

  PrintField("rho");
  PrintField("u");
  PrintField("p");
  PrintField("temperature");
  PrintField("e");
  PrintField("radE");

  FieldsToU(scalar_fields.at("rho"),
            scalar_fields.at("u"),
            scalar_fields.at("v"),
            scalar_fields.at("w"),
            scalar_fields.at("e"),
            U_n);

  rad_E_n = scalar_fields.at("radE");

  const auto num_timesteps = basic_options("max_timesteps").IntegerValue();
  const auto delta_t_max   = basic_options("maximum_dt"   ).FloatValue();
  const auto CFL           = basic_options("CFL"          ).FloatValue();
  const auto max_time      = basic_options("max_time"     ).FloatValue();

  const auto& CCL   = scalar_fields.at("cell_char_length");
  const auto& gamma = scalar_fields.at("gamma");
  const auto& Cv    = scalar_fields.at("Cv");

  auto kappa_s_function = basic_options("kappa_s_function").StringValue();
  auto kappa_a_function = basic_options("kappa_a_function").StringValue();

  const size_t num_local_nodes = grid->local_cells.size();

  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, gamma, CCL, CFL);
    const double dt       = std::min(delta_t_max, dt_hydro);

    chi::log.Log() << "Timestep " << n << " with dt=" << dt << " t=" << time;

    //================================= Compute gradients at n
    //Populates grad_U with slope limiters
    auto grad_U_n     = ComputeUGradients(U_n,*grid,*fv,bc_settings);
    auto grad_rad_E_n = ComputeRadEGradients(rad_E_n,*grid,*fv,bc_settings);

    //=================================== Compute kappas at n
    VecDbl kappa_a_n(num_local_nodes, 0.0);
    VecDbl kappa_t_n(num_local_nodes, 0.0);
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_n,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_n,kappa_t_n);

    Predictor(bc_settings, kappa_a_n, kappa_t_n,
              dt,
              U_n,grad_U_n,U_nph,
              rad_E_n,grad_rad_E_n,rad_E_nph);

    //================================= Compute gradients at nph
    //Populates grad_U with slope limiters
    auto grad_U_nph     = ComputeUGradients(U_nph,*grid,*fv,bc_settings);
    auto grad_rad_E_nph = ComputeRadEGradients(rad_E_nph,*grid,*fv,bc_settings);

    //=================================== Compute kappas at nph
    VecDbl kappa_a_nph(num_local_nodes, 0.0);
    VecDbl kappa_t_nph(num_local_nodes, 0.0);
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_nph,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_nph,kappa_t_nph);

    Corrector(bc_settings, kappa_a_n, kappa_t_n,
              kappa_a_nph, kappa_t_nph,
              dt,
              U_n,U_nph,grad_U_nph,U_np1,
              rad_E_n,rad_E_nph,grad_rad_E_nph,rad_E_np1);

    U_n = U_np1;
    rad_E_n = rad_E_np1;

    time += dt;
    if (max_time > 0.0)
      if (time >= max_time) break;
  }//for n

  UToFields(U_n,
            scalar_fields.at("rho"),
            scalar_fields.at("u"),
            scalar_fields.at("v"),
            scalar_fields.at("w"),
            scalar_fields.at("e"),
            scalar_fields.at("p"),
            scalar_fields.at("gamma"));

  scalar_fields.at("radE") = rad_E_n;

  auto& temperature = scalar_fields.at("temperature");

  const auto& e  = scalar_fields.at("e");

  for (uint64_t c=0; c<num_nodes_local; ++c)
    temperature[c] = e[c] / Cv[c];

  PrintField("rho");
  PrintField("u");
  PrintField("p");
  PrintField("e");
  PrintField("radE");

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}