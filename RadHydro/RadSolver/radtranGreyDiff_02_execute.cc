#include "radtranGreyDiff.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "radtran_structs.h"

void chi_radhydro::RadTranGreyDiffusion::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  std::vector<UVector>&             U_n = U_old;
  std::vector<UVector>              U_nph(num_nodes_local);
  std::vector<UVector>              U_nph_star(num_nodes_local);
  std::vector<UVector>              U_np1(num_nodes_local);

  std::vector<double>&              rad_E_n = rad_E;
  std::vector<double>               rad_E_nph(num_nodes_local);
  std::vector<double>               rad_E_nph_star(num_nodes_local);
  std::vector<double>               rad_E_np1(num_nodes_local);

  FieldsToU(U_n, *grid, rho,u,v,w,e);

  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double dt_hydro = ComputeCourantLimitDelta_t(U_n, gamma,
                                                       cell_char_length, C_cfl);
    const double dt            = std::min(delta_t_max, dt_hydro);

    chi::log.Log() << "Timestep " << n << " with dt=" << dt << " t=" << time;

    //================================= Compute gradients
    //Populates grad_U with slope limiters
    auto grad_U_n     = ComputeUGradients(U_n);
    auto grad_rad_E_n = ComputeRadEGradients(rad_E_n);

    Predictor(dt,
              U_n,
              U_nph,
              grad_U_n,
              rad_E_n,
              rad_E_nph,
              grad_rad_E_n);

    //================================= Compute gradients
    //Populates grad_U with slope limiters
    auto grad_U_nph     = ComputeUGradients(U_n);
    auto grad_rad_E_nph = ComputeRadEGradients(rad_E_n);

    Corrector(dt,
              U_n,
              U_nph,
              U_np1,
              grad_U_nph,
              rad_E_n,
              rad_E_nph,
              rad_E_np1,
              grad_rad_E_nph);

    U_n = U_np1;
    UToFields(U_n,*grid,rho,u,v,w,e,p,gamma);

    time += dt;
    if (t_max >= 0.0 and time >= t_max)
      break;
  }//for n timesteps

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}

