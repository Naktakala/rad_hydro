#ifndef RADHYDRO_LUA_SOLVER_BASE_H
#define RADHYDRO_LUA_SOLVER_BASE_H

#include "ChiLua/chi_lua.h"

#include "../solver_base.h"

namespace chi_radhydro::solver_base_lua_utils
{
  //00
  void RegisterLuaEntities(lua_State* L);

  //01
  int chiRadHydroSolverSetScalarFieldWithLV(lua_State* L);
  //02
  int chiRadHydroSetBCSetting(lua_State* L);
}//chi_radhydro::solver_base_lua_utils

#endif //RADHYDRO_LUA_SOLVER_BASE_H
