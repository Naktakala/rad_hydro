#include "rad_hydro.h"

#include "chi_runtime.h"
#include "chi_log.h"

chi_radhydro::UVector chi_radhydro::
HLLC_RiemannSolve(const UVector &U_L_raw,
                  const UVector &U_R_raw,
                  const double gamma,
                  const Vec3& n_f,
                  bool verbose/*=false*/)
{
  const auto T_Tinv = MakeTandTinv(n_f);

  const auto& T = T_Tinv.first;
  const auto& Tinv = T_Tinv.second;

//  const auto T = chi_radhydro::MakeTransformationMatrix(n_f);
//  const auto Tinv = T.Inverse();

  const UVector U_L = T*U_L_raw;
  const UVector U_R = T*U_R_raw;

  const double p_L = IdealGasPressureFromCellU(U_L, gamma);
  const double p_R = IdealGasPressureFromCellU(U_R, gamma);

  const FVector F_L = MakeFNoTransform(U_L, p_L);
  const FVector F_R = MakeFNoTransform(U_R, p_R);

  const double rho_L = U_L[RHO];
  const double rho_R = U_R[RHO];

  const double u_L = U_L[RHO_U] / rho_L;
  const double u_R = U_R[RHO_U] / rho_R;

  const double a_L = sqrt(gamma*p_L/rho_L);
  const double a_R = sqrt(gamma*p_R/rho_R);

  const double S_L = std::min(u_L-a_L, u_R-a_R);
  const double S_R = std::max(u_L+a_L, u_R+a_R);

  const double S_star = (p_R-p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R))/
                        (rho_L*(S_L-u_L) - rho_R*(S_R-u_R));

  const UVector D_star({0,1,0,0,S_star});

  const double p_star_L = p_L + rho_L*(S_L-u_L)*(S_star-u_L);
  const double p_star_R = p_R + rho_R*(S_R-u_R)*(S_star-u_R);

  const UVector U_star_L = (S_L*U_L - F_L + p_star_L*D_star)/(S_L - S_star);
  const UVector U_star_R = (S_R*U_R - F_R + p_star_R*D_star)/(S_R - S_star);

  const FVector F_star_L = MakeFNoTransform(U_star_L, p_star_L);
  const FVector F_star_R = MakeFNoTransform(U_star_R, p_star_R);

//  const FVector F_star_L =
//    (S_star*(S_L*U_L-F_L) + S_L*(p_L+rho_L*(S_L-u_L)*(S_star-u_L))*D_star)/
//    (S_L-S_star);
//  const FVector F_star_R =
//    (S_star*(S_R*U_R-F_R) + S_R*(p_R+rho_L*(S_R-u_R)*(S_star-u_R))*D_star)/
//    (S_R-S_star);

  bool failed = false;
  FVector F_hllc;
  if      (S_L >= 0)                 F_hllc = F_L;
  else if (S_L <=0 and S_star >= 0)  F_hllc = F_star_L;
  else if (S_star <= 0 and S_R >= 0) F_hllc = F_star_R;
  else if (S_R <= 0)                 F_hllc = F_R;
  else
    failed = true;


  if (verbose or failed)
  {
    chi::log.Log() << "------------";
    chi::log.Log() << "n_f       :" << n_f.PrintStr();

    chi::log.Log() << "u_L       :" << u_L;
    chi::log.Log() << "u_R       :" << u_R;
    chi::log.Log() << "-";
    chi::log.Log() << "p_L       :" << p_L;
    chi::log.Log() << "p_R       :" << p_R;
    chi::log.Log() << "-";
    chi::log.Log() << "a_L       :" << a_L;
    chi::log.Log() << "a_R       :" << a_R;
    chi::log.Log() << "-";
    chi::log.Log() << "S_L       :" << S_L;
    chi::log.Log() << "S_star    :" << S_star;
    chi::log.Log() << "S_R       :" << S_R;
    chi::log.Log() << "-";
    chi::log.Log() << "U_L       :" << PrintU(U_L);
    chi::log.Log() << "U_R       :" << PrintU(U_R);
    chi::log.Log() << "-";
    chi::log.Log() << "F_L       :" << PrintU(F_L     );
    chi::log.Log() << "F_star_L  :" << PrintU(F_star_L);
    chi::log.Log() << "F_star_R  :" << PrintU(F_star_R);
    chi::log.Log() << "F_R       :" << PrintU(F_R     );
    chi::log.Log() << "-";
    chi::log.Log() << "F_hllc    :" << PrintU(F_hllc  );
    chi::log.Log() << "TinvF_hllc:" << PrintU(Tinv * F_hllc  );
  }

  return Tinv * F_hllc;
}