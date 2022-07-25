#include "rad_hydro.h"

#include "ChiMath/chi_math.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Make a 3x3 rotation vector based on an axis and an angle.*/
MatDbl chi_radhydro::MakeRotationMatrix(const Vec3 &axis, double theta)
{
  const auto& ax = axis;
  MatDbl A = {{ 0   ,-ax.z, ax.y},
              { ax.z,    0,-ax.x},
              {-ax.y, ax.x,    0}};
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
/**Makes a 5x5 transformation matrix (for a U-vector)
 * based on a normal.*/
chi_math::MatrixNXxNX<5,double> chi_radhydro::
MakeTransformationMatrix(const Vec3 &n)
{
  constexpr double epsilon = 1.0e-12;
  Vec3 ax(0.0,1.0,0.0);

  const double n_dot_ihat = n.Dot(ihat);
  if (std::fabs(n_dot_ihat)<(1.0-epsilon))
    ax = n.Cross(ihat).Normalized();

  const double theta = acos(n_dot_ihat);

  const MatDbl R = MakeRotationMatrix(ax,theta);
  chi_math::MatrixNXxNX<5,double> T;
  T.SetDiagonalVec({1,0,0,0,1});

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      T.SetIJ(i+1, j+1, R[i][j]);

  return T;
}

chi_radhydro::TMatrixPair chi_radhydro::MakeTandTinv(const Vec3 &n)
{
  const double n_dot_khat = n.Dot(khat);
  if (std::fabs(n_dot_khat) > 0.9999)
  {
    TMatrix T,Tinv;
    if (n_dot_khat > 0)
    {
      T    = {{ 1, 0, 0, 0, 0},
              { 0, 0, 0, 1, 0},
              { 0, 0, 1, 0, 0},
              { 0,-1, 0, 0, 0},
              { 0, 0, 0, 0, 1}};
      Tinv = {{ 1, 0, 0, 0, 0},
              { 0, 0, 0,-1, 0},
              { 0, 0, 1, 0, 0},
              { 0, 1, 0, 0, 0},
              { 0, 0, 0, 0, 1}};
    }
    else
    {
      T    = {{ 1, 0, 0, 0, 0},
              { 0, 0, 0,-1, 0},
              { 0, 0, 1, 0, 0},
              { 0, 1, 0, 0, 0},
              { 0, 0, 0, 0, 1}};
      Tinv = {{ 1, 0, 0, 0, 0},
              { 0, 0, 0, 1, 0},
              { 0, 0, 1, 0, 0},
              { 0,-1, 0, 0, 0},
              { 0, 0, 0, 0, 1}};
    }
    return std::make_pair(T,Tinv);
  }
  else
  {
    auto T    = MakeTransformationMatrix(n);
    auto Tinv = T.Inverse();

    return std::make_pair(T,Tinv);
  }
}

chi_radhydro::FVector chi_radhydro::
  MakeFNoTransform(const UVector &U, double pressure)
{
  const double rho   = U[UVectorEntries::RHO];
  const double rho_u = U[UVectorEntries::RHO_U];

  const double u = rho_u/rho;

  const UVector D({0.0,1.0,0.0,0.0,u});

  return u*U + pressure*D;
}

//###################################################################
/**Makes a flux vector, F, from a U vector.*/
chi_radhydro::FVector chi_radhydro::
MakeF(const UVector& U, double pressure, const Vec3& n_f/*=Vec3(0,0,1)*/)
{
  const auto T_Tinv = MakeTandTinv(n_f);

  const auto& T = T_Tinv.first;
  const auto& Tinv = T_Tinv.second;

//  const auto T = chi_radhydro::MakeTransformationMatrix(n_f);
//  const auto Tinv = T.Inverse();

  const UVector U_t = T * U;

  const double rho   = U_t[0];
  const double rho_u = U_t[1];

  const double u = rho_u/rho;

  const UVector D({0.0,1.0,0.0,0.0,u});

  return Tinv * (u*U_t + pressure*D);
}

//###################################################################
/**Makes a flux vector, F, from a U vector.*/
chi_radhydro::FVector chi_radhydro::
MakeFWithRadE(const UVector& U, const double pressure, const double radE,
              const Vec3& n_f/*=Vec3(0,0,1)*/)
{
  const auto T    = MakeTransformationMatrix(n_f);
  const auto Tinv = T.Inverse();

  const UVector U_t = T * U;

  const double rho   = U_t[0];
  const double rho_u = U_t[1];

  const double u = rho_u/rho;

  const UVector D({0.0,1.0,0.0,0.0,u});
  const UVector e_1_radE({0.0,radE/3.0,0.0,0.0,(4.0/3.0)*radE*u});

  return Tinv * (u*U_t + pressure*D + e_1_radE);
}