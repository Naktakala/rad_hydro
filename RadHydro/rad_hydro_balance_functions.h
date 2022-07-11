#ifndef RADHYDRO_RAD_HYDRO_BALANCE_FUNCTIONS_H
#define RADHYDRO_RAD_HYDRO_BALANCE_FUNCTIONS_H

#include "rad_hydro.h"

#include "ChiMath/chi_math.h"
#include "ChiMath/SerialNewtonIteration/serial_newton_iteration.h"

namespace chi_radhydro
{

  class RadHydroBalanceFunction : public chi_math::NonLinearFunction
  {
  private:
    const double Cv;
    const double gamma;
    const double k1; //=gamma*Cv
    const double k2; //=(gamma-1)*Cv
    const double k3; //=a/3.0
    const VecDbl LHS;

  private:
    VecDbl OneEquationSide(const VecDbl& x) const
    {
      const double& rho = x[0];
      const double& T   = x[1];
      const double& u   = x[2];

      return {rho*u,
              rho*u*u + k2*rho*T + k3*pow(T,4.0),
              0.5*rho*pow(u, 3.0) + k1*rho*T*u + (4.0/3.0)*a*pow(T,4.0)*u};
    }

  public:
    RadHydroBalanceFunction(double in_Cv, double in_gamma, double in_rho0,
                            double in_T0, double in_u0) :
      Cv(in_Cv), gamma(in_gamma),
      k1(gamma*Cv), k2((gamma-1.0)*Cv),
      k3(a/3.0),
      LHS(OneEquationSide({in_rho0, in_T0, in_u0}))
    {}

    VecDbl F(const VecDbl& x) const override
    {
      using namespace chi_math;
      return OneEquationSide(x) - LHS;
    }

    MatDbl J(const VecDbl& x) const override
    {
      const double& rho = x[0];
      const double& T   = x[1];
      const double& u   = x[2];

      MatDbl J =
        {{u             ,                            0, rho},
         {u*u + k2*T + 0, 0 + k2*rho + 4*k3*pow(T,3.0), 2*rho*u + 0 + 0},
           {0.5*pow(u,3.0) + k1*T*u + 0,
              0 + k1*rho*u + 4*(4.0/3.0)*a*pow(T,3.0)*u,
                0.5*3*rho*u*u + k1*rho*T + (4.0/3.0)*a*pow(T,4.0)}};

      return J;
//      const double eps = 1.0e-10;
//
//      auto FD = [&x,this,&eps](size_t i, size_t j)
//      {
//        using namespace chi_math;
//        VecDbl x_plus_eps = x; x_plus_eps[j] += eps;
//
//        auto dF = F(x_plus_eps) - F(x);
//
//        return dF[i]/eps;
//      };
//
//      MatDbl J = {{FD(0,0),FD(0,1),FD(0,2)},
//                  {FD(1,0),FD(1,1),FD(1,2)},
//                  {FD(2,0),FD(2,1),FD(2,2)}};
//
//      return J;
    }
  };

}//namespace chi_radhydro

#endif //RADHYDRO_RAD_HYDRO_BALANCE_FUNCTIONS_H
