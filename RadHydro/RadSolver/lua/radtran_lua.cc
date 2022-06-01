#include "radtran_lua.h"

#include "../radtranDO.h"
#include "chi_runtime.h"

//###################################################################
/**Creates a Discrete Ordinates Radiative Transfer solver.
\param solver_name string (Optional). Text name of the solver. [Default=RadTranDO]
\return solver_handle int
*/
int chi_radtran::lua_utils::chiCreateRadTranDOSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string solver_name = "RadTranDO";
  if (num_args >= 1)
  {
    LuaCheckNilValue(fname, L, 1);
    LuaCheckStringValue(fname, L, 1);

    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<RadTranDO>(solver_name);

  chi::solver_stack.push_back(new_solver);
  const size_t handle = chi::solver_stack.size()-1;
  lua_pushinteger(L, static_cast<lua_Integer>(handle));

  return 1;
}

void chi_radtran::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lua_register(L, "chiCreateRadTranDOSolver", chiCreateRadTranDOSolver);

}