#ifndef RADTRAN_LUA_H
#define RADTRAN_LUA_H

#include "ChiLua/chi_lua.h"

namespace chi_radtran::lua_utils
{
  int chiCreateRadTranDOSolver(lua_State* L);

  void RegisterLuaEntities(lua_State* L);
}//namespace chi_radtran::lua_utils

#endif // RADTRAN_LUA_H
