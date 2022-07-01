#include "rhsolverA.h"

void chi_radhydro::SolverA_GDCN::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

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

  //======================================== Get options
  const auto num_timesteps = basic_options("max_timesteps").IntegerValue();
  const auto delta_t_max   = basic_options("maximum_dt"   ).FloatValue();
  const auto CFL           = basic_options("CFL"          ).FloatValue();
  const auto max_time      = basic_options("max_time"     ).FloatValue();

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
  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, gamma, CCL, CFL);
    const double dt       = std::min(delta_t_max, dt_hydro);

    chi::log.Log() << "Timestep " << n
                   << " with dt=" << dt
                   << " t=" << time
                   << " shock_loc: " << shock_location;

    //================================= Compute gradients at n
    // With slope limiters
    auto grad_U_n     = ComputeUGradients(U_n,*grid,*fv,bc_settings);
    auto grad_rad_E_n = ComputeRadEGradients(rad_E_n,*grid,*fv,bc_settings);

    //================================= Compute kappas at n
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_n,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_n,kappa_t_n);

    //================================= Predictor
    Predictor(bc_settings, kappa_a_n, kappa_t_n,
              dt,
              U_n,grad_U_n,U_nph,
              rad_E_n,grad_rad_E_n,rad_E_nph);

    //================================= Compute gradients at nph
    // With slope limiters
    auto grad_U_nph     = ComputeUGradients(U_nph,*grid,*fv,bc_settings);
    auto grad_rad_E_nph = ComputeRadEGradients(rad_E_nph,*grid,*fv,bc_settings);

    //=================================== Compute kappas at nph
    chi_radhydro::ComputeCellKappas(*grid, Cv, U_nph,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_nph,kappa_t_nph);

    //================================= Corrector
    Corrector(bc_settings, kappa_a_n, kappa_t_n,
              kappa_a_nph, kappa_t_nph,
              dt,
              U_n,U_nph,grad_U_nph,U_np1,
              rad_E_n,rad_E_nph,grad_rad_E_nph,rad_E_np1);

    //================================= Copy new solution to old
    U_n = U_np1;
    rad_E_n = rad_E_np1;

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

    time += dt;
    if (max_time > 0.0)
      if (time >= max_time) break;
  }//for n
  //######################################## Done iterating

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