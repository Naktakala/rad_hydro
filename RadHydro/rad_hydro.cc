#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/chi_math.h"

typedef chi_math::VectorN<5> UVector;
typedef UVector FVector;

#include <iostream>
#include <algorithm>
#include <iomanip>

void chi_radhydro::TestFunction()
{
  std::cout << "Hello from test function\n";
}

//###################################################################
/**Passes the fields-data to a cell U vectors.*/
void chi_radhydro::FieldsToU(std::vector<UVector> &pU,
                             const chi_mesh::MeshContinuum& grid,
                             const std::vector<double>& rho,
                             const std::vector<double>& u,
                             const std::vector<double>& v,
                             const std::vector<double>& w,
                             const std::vector<double>& e)
{
  if (pU.size() != grid.local_cells.size())
    throw std::logic_error("chi_radhydro::FieldsToU used with incompatible"
                           " U-vector dimension.");

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;
    pU[c](0) = rho[c];
    pU[c](1) = rho[c]*u[c];
    pU[c](2) = rho[c]*v[c];
    pU[c](3) = rho[c]*w[c];
    pU[c](4) = rho[c]*(0.5*u[c]*u[c] + 0.5*v[c]*v[c] + 0.5*w[c]*w[c] + e[c]);
  }//for cell
}

//###################################################################
/**Computes the source speed for an ideal gas.*/
double chi_radhydro::
SoundSpeed(double gamma, double p, double rho)
{
  return sqrt(gamma*p/rho);
}

//###################################################################
/**Computes the Courant limit time step.*/
double chi_radhydro::
  ComputeCourantLimitDelta_t(const chi_mesh::MeshContinuum& grid,
                             const std::vector<double>& rho,
                             const std::vector<double>& u,
                             const std::vector<double>& v,
                             const std::vector<double>& w,
                             const std::vector<double>& p,
                             const std::vector<double>& gamma,
                             const std::vector<double>& cell_char_length,
                             const double C_cfl)
{
  double delta_t_hydro = 100.0;
  for (const auto& cell : grid.local_cells)
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

//###################################################################
/**Min mod operator.*/
double chi_radhydro::MinMod(const std::vector<double> &a,
                            bool verbose/*=false*/)
{
  if (verbose)
  {
    std::stringstream outp;
    outp << "MinMod: {";
    for (auto a_i : a) outp << a_i <<",";
    outp << "} ";
    outp << "Signs: {";
    for (auto a_i : a) outp << std::signbit(a_i) <<",";
    outp << "}";

    chi::log.Log() << outp.str();
  }

  if (a.empty()) return 0.0;

  const auto same_sign = std::signbit(a.front());

  for (double a_i : a)
    if (std::signbit(a_i) != same_sign) return 0.0;

  const bool all_negative = bool(same_sign);

  if (not all_negative) return *std::min_element(a.begin(), a.end());
  else                  return *std::max_element(a.begin(), a.end());


}

//###################################################################
/**Applies minmod limiter to a UVector*/
UVector chi_radhydro::MinModU(const std::vector<UVector> &vec_of_U,
                                        bool verbose/*=false*/)
{
  if (vec_of_U.empty())
    throw std::logic_error("chi_hydro::CompInFFlow::MinModU used with no list.");

  UVector minmod_val({0,0,0,0,0});

  const size_t num_elements = vec_of_U.size();
  for (size_t i=0; i<5; ++i)
  {
    std::vector<double> minmod_args(num_elements);

    for (size_t k=0; k<num_elements; ++k)
      minmod_args[k] = vec_of_U[k][static_cast<int>(i)];

    minmod_val(static_cast<int>(i)) = MinMod(minmod_args,verbose);
  }

  return minmod_val;
}

//###################################################################
/**Extrapolates cell-centered U to a face.*/
UVector chi_radhydro::
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
/**Makes a flux vector, F, from a U vector.*/
FVector chi_radhydro::
MakeF(const UVector& U, double pressure, const Vec3& n_f/*=Vec3(1,0,0)*/)
{
  const auto T    = MakeTransformationMatrix(n_f);
  const auto Tinv = T.Inverse();

  const UVector U_t = T * U;

  const double rho   = U_t[0];
  const double rho_u = U_t[1];

  const double u = rho_u/rho;

  const UVector D({0.0,1.0,0.0,0.0,u});

  return Tinv * (u*U_t + pressure*D);
}

//###################################################################
/**Makes a transformation matrix based on a normal.*/
chi_math::MatrixNXxNX<5,double> chi_radhydro::
MakeTransformationMatrix(const Vec3 &n)
{
  constexpr double epsilon = 1.0e-12;
  const Vec3 ihat(1,0,0);
  Vec3 a(0.0,1.0,0.0);

  const double n_dot_ihat = n.Dot(ihat);
  if (std::fabs(n_dot_ihat)<(1.0-epsilon))
    a = n.Cross(ihat).Normalized();

  const double theta = acos(n_dot_ihat);

  const MatDbl R = MakeRotationMatrix(a,theta);
  chi_math::MatrixNXxNX<5,double> T;
  T.SetDiagonalVec({1,0,0,0,1});

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      T.SetIJ(i+1, j+1, R[i][j]);

  return T;
}

//###################################################################
/**Make a rotation vector based on an axis and an angle.*/
MatDbl chi_radhydro::MakeRotationMatrix(const Vec3 &axis, double theta)
{
  const auto& a = axis;
  MatDbl A = {{ 0  ,-a.z, a.y},
              { a.z,   0,-a.x},
              {-a.y, a.x,   0}};
  MatDbl I = {{1,0,0},
              {0,1,0},
              {0,0,1}};
  const auto AA = chi_math::MatMul(A,A);
  const auto C0 = chi_math::MatMul(AA, 1.0-cos(theta));
  const auto C1 = chi_math::MatMul(A, sin(theta));
  const auto C2 = chi_math::MatAdd(C0, C1);

  return chi_math::MatAdd(I, C2);
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

UVector chi_radhydro::
HLLC_RiemannSolve(const UVector &U_L_raw,
                  const UVector &U_R_raw,
                  const double p_L,
                  const double p_R,
                  const double gamma_L,
                  const double gamma_R,
                  const Vec3& n_f,
                  bool verbose/*=false*/)
{
  const auto T    = MakeTransformationMatrix(n_f);
  const auto Tinv = T.Inverse();

  const UVector U_L = T*U_L_raw;
  const UVector U_R = T*U_R_raw;

  const FVector F_L = MakeF(U_L, p_L);
  const FVector F_R = MakeF(U_R, p_R);

  const double rho_L = U_L[0];
  const double rho_R = U_R[0];

  const double u_L = U_L[1]/rho_L;
  const double u_R = U_R[1]/rho_R;

  const double a_L = sqrt(gamma_L*p_L/rho_L);
  const double a_R = sqrt(gamma_R*p_R/rho_R);

  const double S_L = std::min(u_L-a_L, u_R-a_R);
  const double S_R = std::max(u_L+a_L, u_R+a_R);

  const double S_star = (p_R-p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R))/
                        (rho_L*(S_L-u_L) - rho_R*(S_R-u_R));

  const UVector D_star({0,1,0,0,S_star});

  const FVector F_star_L =
    (S_star*(S_L*U_L-F_L) + S_L*(p_L+rho_L*(S_L-u_L)*(S_star-u_L))*D_star)/
    (S_L-S_star);
  const FVector F_star_R =
    (S_star*(S_R*U_R-F_R) + S_R*(p_R+rho_L*(S_R-u_R)*(S_star-u_R))*D_star)/
    (S_R-S_star);

  FVector F_hllc;
  if (S_L >= 0)                 F_hllc = F_L;
  if (S_L <=0 and S_star >= 0)  F_hllc = F_star_L;
  if (S_star <= 0 and S_R >= 0) F_hllc = F_star_R;
  if (S_R <= 0)                 F_hllc = F_R;

  if (verbose)
  {
    chi::log.Log() << "S_L   :" << S_L;
    chi::log.Log() << "S_star:" << S_star;
    chi::log.Log() << "S_R   :" << S_R;

    chi::log.Log() << "U_L:" << PrintU(U_L);
    chi::log.Log() << "U_R:" << PrintU(U_R);

    chi::log.Log() << "F_L     :" << PrintU(F_L     );
    chi::log.Log() << "F_star_L:" << PrintU(F_star_L);
    chi::log.Log() << "F_star_R:" << PrintU(F_star_R);
    chi::log.Log() << "F_R     :" << PrintU(F_R     );
  }

  return Tinv * F_hllc;
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
      output << std::setprecision(3) << std::scientific << val << " ";
    else
      output << std::setprecision(3) << std::scientific << 0.0 << " ";
  }
  output << "]";

  return output.str();
}

//###################################################################
/**Passes cell U vectors to fields data.*/
void chi_radhydro::
UToFields(const std::vector<UVector> &pU,
          const chi_mesh::MeshContinuum& grid,
          std::vector<double>& rho,
          std::vector<double>& u,
          std::vector<double>& v,
          std::vector<double>& w,
          std::vector<double>& e,
          std::vector<double>& p,
          const std::vector<double>& gamma)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;

    rho[c] = pU[c][0];
    u[c]   = pU[c][1]/rho[c];
    v[c]   = pU[c][2]/rho[c];
    w[c]   = pU[c][3]/rho[c];

    const double E = pU[c][4];
    e[c]   = E/rho[c] - 0.5*u[c]*u[c] - 0.5*v[c]*v[c] - 0.5*w[c]*w[c];
    p[c]   = (gamma[c]-1.0)*rho[c]*e[c];
  }//for cell
}