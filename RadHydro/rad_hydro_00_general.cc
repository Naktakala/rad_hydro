#include "rad_hydro.h"

//###################################################################
/**Very convenient to get and reset a timer.*/
double chi_radhydro::GetResetTimer(chi_objects::ChiTimer &timer)
{
  double time = timer.GetTime();
  timer.Reset();
  return time;
}

//###################################################################
/**Computes the source speed for an ideal gas.*/
double chi_radhydro::
SoundSpeed(double gamma, double p, double rho)
{
  return sqrt(gamma*p/rho);
}

//###################################################################
/**Computes the courant limited timestep size for the given
 * hydrodynamic state U.*/
double chi_radhydro::
  ComputeCourantLimitDelta_t(const std::vector<UVector> &vecU,
                             const double gamma,
                             const std::vector<double> &cell_char_length,
                             double C_cfl)
{

  const size_t N = vecU.size();
  double delta_t_hydro = 100.0;
  for (size_t c=0; c<N; ++c)
  {
    const auto& U_c = vecU[c];

    const double u_norm_c = chi_radhydro::VelocityFromCellU(U_c).Norm();
    const double p_c      = chi_radhydro::IdealGasPressureFromCellU(U_c, gamma);
    const double rho_c    = U_c[0];
    const double a_c      = SoundSpeed(gamma,p_c,rho_c);
    const double L_c      = cell_char_length[c];

    const double u_norm = u_norm_c + a_c;
    if (c == 0)
      delta_t_hydro =          2.0 * C_cfl * L_c / u_norm;
    else
      delta_t_hydro = std::min(2.0 * C_cfl * L_c / u_norm, delta_t_hydro);
  }//for cell

  return delta_t_hydro;
}

//###################################################################
/**Makes a Uvector from a boundary condition.*/
chi_radhydro::UVector chi_radhydro::MakeUFromBC(const BCSetting& bc_setting,
                                                const UVector& U_default)
{
  if (bc_setting.type == BCType::TRANSMISSIVE) return U_default;

  UVector U;
  double rho = bc_setting.values[0];
  double u   = bc_setting.values[1];
  double v   = bc_setting.values[2];
  double w   = bc_setting.values[3];
  double e   = bc_setting.values[4];

  Vec3 vel = Vec3(u,v,w);

  double E = rho*(0.5*vel.NormSquare() + e);

  U(0) = rho;
  U(1) = rho*u;
  U(2) = rho*v;
  U(3) = rho*w;
  U(4) = E;

  return U;
}

//###################################################################
/**Makes radiation energy from a boundary condition.*/
double chi_radhydro::MakeRadEFromBC(const BCSetting& bc_setting,
                                    const double radE_default)
{
  if (bc_setting.type == BCType::TRANSMISSIVE) return radE_default;

  return bc_setting.values[6];
}