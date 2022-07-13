#include "rad_hydro.h"

#include <iomanip>

//###################################################################
/**Extrapolates cell-centered U to a face.*/
chi_radhydro::UVector chi_radhydro::
UplusDXGradU(const UVector &U,
             const Vec3 &dx,
             const GradUTensor &grad_U)
{
  UVector U_f = U + dx.x * grad_U[0] +
                dx.y * grad_U[1] +
                dx.z * grad_U[2];

  return U_f;
}

//###################################################################
/**Computes the cell pressure from the ideal gas law for a single cell.*/
double chi_radhydro::
IdealGasPressureFromCellU(const UVector &U, double gamma)
{
  const double rho_c = U[0];
  const double u_c   = U[1]/rho_c;
  const double v_c   = U[2]/rho_c;
  const double w_c   = U[3]/rho_c;

  const double E   = U[4];
  const double e_c = E/rho_c - 0.5*u_c*u_c - 0.5*v_c*v_c - 0.5*w_c*w_c;
  const double p_c = (gamma-1.0)*rho_c*e_c;

  return p_c;
}

//###################################################################
/**Computes the temperature from the ideal gas law for a single cell.*/
double chi_radhydro::
IdealGasTemperatureFromCellU(const UVector &U, double C_v)
{
  const double rho_c = U[0];
  const double u_c   = U[1]/rho_c;
  const double v_c   = U[2]/rho_c;
  const double w_c   = U[3]/rho_c;

  const double E   = U[4];
  const double e_c = E/rho_c - 0.5*u_c*u_c - 0.5*v_c*v_c - 0.5*w_c*w_c;

  return e_c/C_v;
}

//###################################################################
/**Computes the internal energy from cell U.*/
double chi_radhydro::
InternalEnergyFromCellU(const UVector &U)
{
  const double rho_c = U[0];
  const double u_c   = U[1]/rho_c;
  const double v_c   = U[2]/rho_c;
  const double w_c   = U[3]/rho_c;

  const double E   = U[4];
  const double e_c = E/rho_c - 0.5*u_c*u_c - 0.5*v_c*v_c - 0.5*w_c*w_c;

  return e_c;
}

//###################################################################
/**Computes the velocity vector from cell U.*/
chi_mesh::Vector3 chi_radhydro::
VelocityFromCellU(const UVector &U)
{
  const double rho_c = U[0];
  const double u_c   = U[1]/rho_c;
  const double v_c   = U[2]/rho_c;
  const double w_c   = U[3]/rho_c;

  return chi_mesh::Vector3(u_c,v_c,w_c);
}

//###################################################################
/**Provides a string representation of a U-vector.*/
std::string chi_radhydro::PrintU(const UVector &U,
                                 const double epsilon/*=1.0e-12*/)
{
  std::stringstream output;
  output << "[";
  for (double val : U.elements)
  {
    if (std::fabs(val) > epsilon)
      output << std::setprecision(8) << std::scientific << val << " ";
    else
      output << std::setprecision(8) << std::scientific << 0.0 << " ";
  }
  output << "]";

  return output.str();
}