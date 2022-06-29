#include "rad_hydro.h"

#include "chi_runtime.h"
#include "ChiConsole/chi_console.h"

//###################################################################
/**This function calls a lua function with two arguments. The first
 * argument is a scalar, and the second argument is an integer.*/
double chi_radhydro::ComputeKappaFromLua(double T, int mat_id, const std::string& lua_fname)
{
  if (lua_fname.empty()) return 0.0;

  double ret_val = 0.0;

  auto& L = chi::console.consoleState;

  int result = lua_getglobal(L,lua_fname.c_str());
  if (result != LUA_TFUNCTION)
  {
    lua_pop(L,1);
    throw std::logic_error(std::string(__FUNCTION__) + ": Tasked with "
    "calling lua-function \"" + lua_fname + "\" but the global definition was not"
    "found.");
  }
  lua_pushnumber(L,T);
  lua_pushinteger(L,mat_id);

  //2 arguments, 1 result, 0=original error object
  if (lua_pcall(L,2,1,0) == 0)
  {
    ret_val = lua_tonumber(L,-1);
  }
  lua_pop(L,1);

  return ret_val;
}