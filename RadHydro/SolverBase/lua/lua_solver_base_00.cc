#include "lua_solver_base.h"

#define REGISTER(x) \
lua_register(L, #x, x)

namespace chi_radhydro::solver_base_lua_utils
{

void RegisterLuaEntities(lua_State* L)
{
  REGISTER(chiRadHydroSolverSetScalarFieldWithLV);
  REGISTER(chiRadHydroSetBCSetting);
}

}//namespace chi_radhydro::solver_base_lua_utils