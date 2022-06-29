#include "rad_hydro.h"

#include "ChiMath/chi_math.h"

//###################################################################
/**Make a 3x3 rotation vector based on an axis and an angle.*/
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
/**Makes a 5x5 transformation matrix (for a U-vector)
 * based on a normal.*/
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
/**Makes a flux vector, F, from a U vector.*/
chi_radhydro::FVector chi_radhydro::
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