#include "rhsolverB.h"

namespace chi_radhydro
{

void SolverB_GDTBDF::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  //======================================== Allocate work data
  const size_t num_local_nodes = grid->local_cells.size();
  std::vector<UVector>  U_n(num_local_nodes);
  std::vector<UVector>  U_npq(num_local_nodes);
  std::vector<UVector>  U_nph(num_local_nodes);
  std::vector<UVector>  U_np3q(num_local_nodes);
  std::vector<UVector>  U_np1(num_local_nodes);

  VecDbl                rad_E_n(num_local_nodes);
  VecDbl                rad_E_npq(num_local_nodes);
  VecDbl                rad_E_nph(num_local_nodes);
  VecDbl                rad_E_np3q(num_local_nodes);
  VecDbl                rad_E_np1(num_local_nodes);

  VecDbl                kappa_a_n(num_local_nodes, 0.0);
  VecDbl                kappa_t_n(num_local_nodes, 0.0);

  VecDbl                kappa_a_npq(num_local_nodes, 0.0);
  VecDbl                kappa_t_npq(num_local_nodes, 0.0);

  VecDbl                kappa_a_nph(num_local_nodes, 0.0);
  VecDbl                kappa_t_nph(num_local_nodes, 0.0);

  VecDbl                kappa_a_np3q(num_local_nodes, 0.0);
  VecDbl                kappa_t_np3q(num_local_nodes, 0.0);

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

  auto kappa_s_function = basic_options("kappa_s_function").StringValue();
  auto kappa_a_function = basic_options("kappa_a_function").StringValue();

  //======================================== Init shock location tracker
  double shock_ref_speed = 0.0;
  for (double uval : scalar_fields.at("w"))
    shock_ref_speed = std::max(shock_ref_speed, uval);

  double shock_location = 0.5;

  //######################################## Start iterations
  double time = 0.0;
  double max_change = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, m_gamma, CCL, CFL);
    const double dt       = std::min(delta_t_max, dt_hydro);

    chi::log.Log() << "Timestep " << n
                   << " with dt=" << dt
                   << " t=" << time
                   << " shock_loc: " << shock_location
                   << " l_inf: " << max_change;

    //########################################### Phase A
    //================================= Compute gradients at n
    // With slope limiters
    auto grad_U_n     = ComputeUGradients(U_n,*grid,*fv,bc_settings);
    auto grad_rad_E_n = ComputeRadEGradients(rad_E_n,*grid,*fv,bc_settings);

    //================================= Compute kappas at n
    chi_radhydro::ComputeCellKappas(*grid, m_Cv, U_n,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_n,kappa_t_n);

    //================================= Predictor
    Predictor(bc_settings, kappa_a_n, kappa_t_n,
              dt/2.0,
              U_n,grad_U_n,U_npq,
              rad_E_n,grad_rad_E_n,rad_E_npq);

    //================================= Compute gradients at nph
    // With slope limiters
    auto grad_U_npq     = ComputeUGradients(U_npq,*grid,*fv,bc_settings);
    auto grad_rad_E_npq = ComputeRadEGradients(rad_E_npq,*grid,*fv,bc_settings);

    //=================================== Compute kappas at nph
    chi_radhydro::ComputeCellKappas(*grid, m_Cv, U_nph,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_npq,kappa_t_npq);

    //================================= Corrector
    Corrector(bc_settings,
              kappa_a_n, kappa_t_n,
              kappa_a_npq, kappa_t_npq,
              dt/2.0,
              U_n, U_npq, grad_U_npq, U_nph,
              rad_E_n, rad_E_npq, grad_rad_E_npq, rad_E_nph);

    //########################################### Phase B
    //================================= Compute gradients at n
    // With slope limiters
    auto grad_U_nph     = ComputeUGradients(U_nph,*grid,*fv,bc_settings);
    auto grad_rad_E_nph = ComputeRadEGradients(rad_E_nph,*grid,*fv,bc_settings);

    //================================= Compute kappas at n
    chi_radhydro::ComputeCellKappas(*grid, m_Cv, U_nph,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_nph,kappa_t_nph);

    //================================= Predictor
    Predictor(bc_settings,
              kappa_a_nph, kappa_t_nph,
              dt/2.0,
              U_nph,grad_U_nph,U_np3q,
              rad_E_nph,grad_rad_E_nph,rad_E_np3q);

    //================================= Compute gradients at nph
    // With slope limiters
    auto grad_U_np3q     = ComputeUGradients(U_np3q,*grid,*fv,bc_settings);
    auto grad_rad_E_np3q  = ComputeRadEGradients(rad_E_np3q,*grid,*fv,bc_settings);

    //=================================== Compute kappas at nph
    chi_radhydro::ComputeCellKappas(*grid, m_Cv, U_np3q,
                                    kappa_s_function, kappa_a_function,
                                    kappa_a_np3q,kappa_t_np3q);

    //================================= Corrector
    CorrectorB(bc_settings,
               kappa_a_n, kappa_t_n,
               kappa_a_nph, kappa_t_nph,
               kappa_a_np3q, kappa_t_np3q,
               dt / 2.0,
               U_n, U_nph, U_np3q, grad_U_np3q, U_np1,
               rad_E_n, rad_E_nph, rad_E_np3q, grad_rad_E_np3q, rad_E_np1);




    //================================= Checking for steady state
    max_change = 0.0;
    for (size_t c=0; c<num_local_nodes; ++c)
    {
      const double radT_c_n   = pow((rad_E_n  [c]/a),0.25);
      const double radT_c_np1 = pow((rad_E_np1[c]/a),0.25);
      const double nodal_change = (radT_c_np1 - radT_c_n)/radT_c_n;

      max_change = std::max(max_change,std::fabs(nodal_change));
    }//for k

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
            m_gamma);

  scalar_fields.at("radE") = rad_E_n;

  //======================================== Compute material temperature
  auto& temperature = scalar_fields.at("temperature");

  const auto& e  = scalar_fields.at("e");

  for (uint64_t c=0; c < num_local_nodes; ++c)
    temperature[c] = e[c] / m_Cv;

  //======================================== Print outputs
  PrintRawOutput("ZRawOutput.txt");

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}

}//namespace chi_radhydro