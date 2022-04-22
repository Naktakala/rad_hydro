#include "compinfflow.h"

//###################################################################
/**Computes the Courant limit time step.*/
double chi_hydro::CompInFFlow::ComputeCourantLimitDelta_t() const
{
  double delta_t_hydro = 100.0;
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t c        = cell.local_id;
    const double   u_norm_c = Vec3(u[c],v[c],w[c]).Norm();
    const double   a_c      = SoundSpeed(gamma[c],p[c],rho[c]);
    const double   L_c      = cell_char_length[c];

    const double u_norm = u_norm_c + a_c;
    if (cell.local_id == 0)
      delta_t_hydro =          2.0 * C_cfl * L_c / u_norm;
    else
      delta_t_hydro = std::min(2.0 * C_cfl * L_c / u_norm, delta_t_hydro);
  }//for cell

  return delta_t_hydro;
}