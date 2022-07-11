#include "rad_hydro.h"

#include "rad_hydro_balance_functions.h"

namespace chi_radhydro
{


std::vector<double> MakePostShockConditionsRH(const double Cv,
                                              const double gamma,
                                              const double rho0,
                                              const double T0,
                                              const double u0,
                                              const double rho1_guess,
                                              const double T1_guess,
                                              const double u1_guess)
{
  RadHydroBalanceFunction rad_hydro_balance_function(Cv,gamma,rho0,T0,u0);

  VecDbl x1_guess = {rho1_guess, T1_guess, u1_guess};

  VecDbl x1 = chi_math::NewtonIteration(
    rad_hydro_balance_function, //Non-linear function
    x1_guess,                   //Initial guess
    100,                        //Max iterations
    1.0e-20,                    //Convergence criteria
    true);                      //Verbosity

  return x1;
}

}//namespace chi_radhydro