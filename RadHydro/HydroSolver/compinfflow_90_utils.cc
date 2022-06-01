#include "compinfflow.h"

#include "chi_log.h"

typedef chi_math::VectorN<5> UVector;
typedef UVector FVector;

//###################################################################
/**Computes the cell pressure from the ideal gas law for a single cell.*/
double chi_hydro::CompInFFlow::
  PressureFromCellU(const UVector &U, double gamma)
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
/**Makes a flux vector, F, from a U vector.*/
FVector chi_hydro::CompInFFlow::
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
/**Computes the source speed for an ideal gas.*/
double chi_hydro::CompInFFlow::
  SoundSpeed(double gamma, double p, double rho)
{
  return sqrt(gamma*p/rho);
}

//###################################################################
/**Passes the fields-data to a cell U vectors.*/
void chi_hydro::CompInFFlow::FieldsToU(std::vector<UVector> &pU)
{
  if (pU.size() != grid->local_cells.size())
    throw std::logic_error("CompInFFlow::FieldsToU used with incompatible"
                           " U-vector dimension.");

  for (const auto& cell : grid->local_cells)
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
/**Passes cell U vectors to fields data.*/
void chi_hydro::CompInFFlow::UToFields(const std::vector<UVector> &pU)
{
  for (const auto& cell : grid->local_cells)
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

//###################################################################
/**Make a rotation vector based on an axis and an angle.*/
MatDbl chi_hydro::CompInFFlow::
  MakeRotationMatrix(const Vec3 &axis, double theta)
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
/**Makes a transformation matrix based on a normal.*/
chi_math::MatrixNXxNX<5,double> chi_hydro::CompInFFlow::
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
/**Min mod operator.*/
double chi_hydro::CompInFFlow::MinMod(const std::vector<double> &a,
                                      bool verbose/*=false*/)
{
  if (verbose)
  {
    std::stringstream outp;
    outp << "MinMod: {";
    for (auto a_i : a) outp << a_i <<",";
    outp << "} ";
    outp << "Signs: {";
    for (auto a_i : a) outp << signbit(a_i) <<",";
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
UVector chi_hydro::CompInFFlow::MinModU(const std::vector<UVector> &vec_of_U,
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
UVector chi_hydro::CompInFFlow::
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
/**Provides a string representation of a U-vector.*/
std::string chi_hydro::CompInFFlow::PrintU(const UVector &U,
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

