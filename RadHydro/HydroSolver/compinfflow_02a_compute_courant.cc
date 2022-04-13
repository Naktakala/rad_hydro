#include "compinfflow.h"

//###################################################################
/**Computes the Courant limit time step.*/
double chi_hydro::CompInFFlow::ComputeCourantLimitDelta_t() const
{
  double delta_t_hydro = 100.0;
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t c = cell.local_id;
    const double u_abs = sqrt(pow(u[c],2.0)+pow(v[c],2.0)+pow(w[c],2.0));
    const double a = SoundSpeed(gamma[c],p[c],rho[c]);
    const double dx = cell_char_length[c];

    const double S_max = u_abs + a;
    delta_t_hydro = std::min(2.0*C_cfl*dx/S_max, delta_t_hydro);
  }//for cell

  return delta_t_hydro;
}