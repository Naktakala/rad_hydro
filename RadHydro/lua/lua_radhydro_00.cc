#include "lua_radhydro.h"

#include "../rad_hydro.h"

namespace chi_radhydro::lua_utils
{

int chiRadHydroMakePostShockConditionsRH(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 8)
    LuaPostArgAmountError(fname, 8, num_args);

  const double Cv         = lua_tonumber(L, 1);
  const double gamma      = lua_tonumber(L, 2);
  const double rho0       = lua_tonumber(L, 3);
  const double T0         = lua_tonumber(L, 4);
  const double u0         = lua_tonumber(L, 5);
  const double rho1_guess = lua_tonumber(L, 6);
  const double T1_guess   = lua_tonumber(L, 7);
  const double u1_guess   = lua_tonumber(L, 8);

  VecDbl x1 = MakePostShockConditionsRH(Cv, gamma, rho0, T0, u0,
                                        rho1_guess, T1_guess, u1_guess);

  lua_pushnumber(L, x1[0]); //rho1
  lua_pushnumber(L, x1[1]); //T1
  lua_pushnumber(L, x1[2]); //u1

  return 3;
}

}//namespace chi_radhydro::lua_utils