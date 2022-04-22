#include "compinfflow_lua.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"

//###################################################################
chi_hydro::CompInFFlow* chi_hydro::compinfflow_lua_utils::
  GetSolverByHandle(int handle, const std::string &calling_function_name)
{
  chi_hydro::CompInFFlow* compinfflow_solver;
  try{
    compinfflow_solver = dynamic_cast<chi_hydro::CompInFFlow*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not compinfflow_solver)
      throw std::logic_error(calling_function_name +
       ": Invalid solver at given handle (" +
       std::to_string(handle) + "). "
       "The solver is not of type chi_hydro::CompInFFlow.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }//catch

  return compinfflow_solver;
}

//###################################################################
/**Registers lua-callable entities.*/
void chi_hydro::compinfflow_lua_utils::RegisterLuaEntities(lua_State *L)
{
  lua_register(L, "chiCreateCompInFFlowSolver", chiCreateCompInFFlowSolver);
  lua_register(L, "chiCompInFFlowSetFieldInitialValue", chiCompInFFlowSetFieldInitialValue);
}

//###################################################################
/**Creates a CompInFFlow solver.

\param solver_name string Text name of the solver.

*/
int chi_hydro::compinfflow_lua_utils::chiCreateCompInFFlowSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  auto& chi_log = ChiLog::GetInstance();

  std::string solver_name = "CompInFFlow";
  if (num_args>=1)
  {
    LuaCheckNilValue(fname, L, 1);
    LuaCheckStringValue(fname, L, 1);

    solver_name = lua_tostring(L, 1);
  }

  chi_log.Log() << "Creating CompInFFlow solver with name \""
                << solver_name << "\".";

  auto solver = new chi_hydro::CompInFFlow(solver_name);

  chi_physics_handler.solver_stack.push_back(solver);

  auto n = static_cast<lua_Integer>(chi_physics_handler.solver_stack.size() - 1);
  lua_pushinteger(L, n);
  return 1;
}
