#include "compinfflow.h"

typedef chi_math::VectorN<5> UVector;

#include "chi_log.h"
extern ChiLog& chi_log;

UVector chi_hydro::CompInFFlow::
  HLLC_RiemannSolve(const UVector &U_L,
                    const UVector &U_R,
                    const double p_L,
                    const double p_R,
                    const double gamma_L,
                    const double gamma_R,
                    bool verbose/*=false*/)
{
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

  const double C0_L = rho_L*(S_L-u_L)/(S_L-S_star);
  const double C0_R = rho_R*(S_R-u_R)/(S_R-S_star);

  const UVector U_star_L({C0_L * 1,
                          C0_L * S_star,
                          C0_L * U_L[2]/rho_L,
                          C0_L * U_L[3]/rho_L,
                          C0_L * ((U_L[4]/rho_L) + (S_star-u_L) * (S_star+p_L/(rho_L*(S_L-u_L))))});

  const UVector U_star_R({C0_R * 1,
                          C0_R * S_star,
                          C0_R * U_R[2]/rho_R,
                          C0_R * U_R[3]/rho_R,
                          C0_R * ((U_R[4]/rho_R) + (S_star-u_R) * (S_star+p_R/(rho_R*(S_R-u_R))))});

  const FVector F_star_L = F_L + S_L * (U_star_L - U_L);
  const FVector F_star_R = F_R + S_R * (U_star_R - U_R);

//  const UVector D_star({0,1,0,0,S_star});

//  const FVector F_star_L =
//    (S_star*(S_L*U_L-F_L) + S_L*(p_L+rho_L*(S_L-u_L)*(S_star-u_L))*D_star)/
//    (S_L-S_star);
//  const FVector F_star_R =
//    (S_star*(S_R*U_R-F_R) + S_R*(p_R+rho_L*(S_R-u_R)*(S_star-u_R))*D_star)/
//    (S_R-S_star);

  FVector F_hllc;
  if (S_L >= 0)                 F_hllc = F_L;
  if (S_L <=0 and S_star >= 0)  F_hllc = F_star_L;
  if (S_star <= 0 and S_R >= 0) F_hllc = F_star_R;
  if (S_R <= 0)                 F_hllc = F_R;

  if (verbose)
  {
    chi_log.Log() << "S_L   :" << S_L;
    chi_log.Log() << "S_star:" << S_star;
    chi_log.Log() << "S_R   :" << S_R;

    chi_log.Log() << "U_Lstar:" << PrintU(U_star_L);
    chi_log.Log() << "U_Rstar:" << PrintU(U_star_R);

    chi_log.Log() << "U_L:" << PrintU(U_L);
    chi_log.Log() << "U_R:" << PrintU(U_R);

    chi_log.Log() << "F_L     :" << PrintU(F_L     );
    chi_log.Log() << "F_star_L:" << PrintU(F_star_L);
    chi_log.Log() << "F_star_R:" << PrintU(F_star_R);
    chi_log.Log() << "F_R     :" << PrintU(F_R     );
  }

  return F_hllc;
}