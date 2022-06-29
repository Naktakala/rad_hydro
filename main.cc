#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "RadHydro/HydroSolver/lua/compinfflow_lua.h"
#include "RadHydro/RadSolver/lua/radtran_lua.h"
#include "RadHydro/SolverBase/lua/lua_solver_base.h"
#include "RadHydro/SolverA_GDCN/lua/lua_rhsolverA.h"

int main(int argc, char* argv[])
{
  chi::log.Log() << "chiRadHydro - Execution started";

  chi::Initialize(argc,argv);

  auto& L = chi::console.consoleState;
  chi_hydro::compinfflow_lua_utils::RegisterLuaEntities(L);
  chi_radhydro::lua_utils::RegisterLuaEntities(L);
  chi_radhydro::solver_base_lua_utils::RegisterLuaEntities(L);

  chi::console.RegisterFunction(
    "chiCreateSolverA", chi_radhydro::solverA_lua_utils::chiCreateSolverA);

  chi::RunBatch(argc,argv);

  chi::Finalize();

  chi::log.Log() << "chiRadHydro - Execution finished";

  return 0;
}