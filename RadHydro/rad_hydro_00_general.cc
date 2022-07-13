#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Computes the source speed for an ideal gas.*/
double chi_radhydro::
SoundSpeed(double gamma, double p, double rho)
{
  return sqrt(gamma*p/rho);
}

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

void chi_radhydro::
  ComputeCellKappas(const chi_mesh::MeshContinuum &grid,
                    const double               Cv,
                    const std::vector<UVector> &U,
                    const std::string &kappa_s_function,
                    const std::string &kappa_a_function,
                    std::vector<double> &kappa_a,
                    std::vector<double> &kappa_t)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;
    const int mat_id = cell.material_id;

    const double T_c = IdealGasTemperatureFromCellU(U[c], Cv);

    const double kappa_s_c = ComputeKappaFromLua(T_c, mat_id, kappa_s_function);
    const double kappa_a_c = ComputeKappaFromLua(T_c, mat_id, kappa_a_function);

    kappa_a[c] = kappa_a_c;
    kappa_t[c] = kappa_s_c + kappa_a_c;
  }//for cell
}

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

double chi_radhydro::MakeRadEFromBC(const BCSetting& bc_setting,
                                    const double radE_default)
{
  if (bc_setting.type == BCType::TRANSMISSIVE) return radE_default;

  return bc_setting.values[6];
}

double chi_radhydro::GetResetTimer(chi_objects::ChiTimer &timer)
{
  double time = timer.GetTime();
  timer.Reset();
  return time;
}

/**Computes Emission-absorption source.
 * Sea = sigma_a c (aT^4 - radE)*/
double chi_radhydro::MakeEmAbsSource(const double sigma_a,
                                     const double T,
                                     const double radE)
{
  return sigma_a * speed_of_light_cmpsh *
         (a * pow(T,4.0) - radE);
}

double chi_radhydro::Make3rdGradRadE(const chi_mesh::Cell &cell,
                                     double V_c,
                                     const std::vector<double> &face_areas,
                                     double radE,
                                     const Vec3 &grad_radE,
                                     const UVector &U,
                                     const std::vector<UVector> &grad_U)
{
  double third_grad_rad_E_dot_u = 0.0;
  const auto& x_c = cell.centroid;

//  const Vec3    u_f     = VelocityFromCellU(U);

  size_t f = 0;
  for (const auto& face : cell.faces)
  {
    const auto& x_f  = face.centroid;
    const auto  x_fc = x_f - x_c;
    const Vec3 A_f = face_areas[f] * face.normal;

    const double  rad_E_f = radE + (x_f - x_c).Dot(grad_radE);
    const UVector U_f     = UplusDXGradU(U, x_fc, grad_U);
    const Vec3    u_f     = VelocityFromCellU(U_f);

    third_grad_rad_E_dot_u += (1/V_c)*(1.0/3) * A_f.Dot(rad_E_f * u_f);
    ++f;
  }//for f
  return third_grad_rad_E_dot_u;
}
