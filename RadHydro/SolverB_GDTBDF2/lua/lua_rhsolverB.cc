#include "lua_rhsolverB.h"

#include "../rhsolverB.h"

#include "chi_runtime.h"
#include "chi_log.h"

int chi_radhydro::solverB_lua_utils::chiCreateSolverB(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string solver_name = "SolverB";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto solver = std::make_shared<SolverB_GDTBDF>(solver_name);

  chi::solver_stack.push_back(solver);
  size_t handle = chi::solver_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(handle));

  chi::log.Log() << solver_name << " created.";

  return 1;
}