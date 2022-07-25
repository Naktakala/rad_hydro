#ifndef RADHYDRO_LUA_RADHYDRO_H
#define RADHYDRO_LUA_RADHYDRO_H

#include "ChiLua/chi_lua.h"

namespace chi_radhydro::lua_utils
{

int chiRadHydroMakePostShockConditionsRH(lua_State* L);
int chiRadHydroMakePostShockConditionsHydroOnly(lua_State* L);
int chiRadHydroTestMFEM(lua_State* L);

}//namespace chi_radhydro

#endif //RADHYDRO_LUA_RADHYDRO_H
