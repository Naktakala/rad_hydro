#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "RadHydro/lua/lua_radhydro.h"
#include "RadHydro/HydroSolver/lua/compinfflow_lua.h"
#include "RadHydro/SolverBase/lua/lua_solver_base.h"
#include "RadHydro/SolverA_GDCN/lua/lua_rhsolverA.h"
#include "RadHydro/SolverB_GDCN_MFEM/lua/lua_rhsolverB.h"
#include "RadHydro/SolverC_SNCN_MFEM/lua/lua_rhsolverC.h"

int main(int argc, char* argv[])
{
  chi::log.Log() << "chiRadHydro - Execution started";

  chi::Initialize(argc,argv);

  auto& L = chi::console.consoleState;
  chi_hydro::compinfflow_lua_utils::RegisterLuaEntities(L);
  chi_radhydro::solver_base_lua_utils::RegisterLuaEntities(L);

  chi::console.RegisterFunction(
    "chiCreateSolverA", chi_radhydro::solverA_lua_utils::chiCreateSolverA);
  chi::console.RegisterFunction(
    "chiCreateSolverB", chi_radhydro::solverB_lua_utils::chiCreateSolverB);
  chi::console.RegisterFunction(
    "chiCreateSolverC", chi_radhydro::solverC_lua_utils::chiCreateSolverC);

  chi::console.RegisterFunction(
    "chiRadHydroMakePostShockConditionsRH",
    chi_radhydro::lua_utils::chiRadHydroMakePostShockConditionsRH);
  chi::console.RegisterFunction(
    "chiRadHydroMakePostShockConditionsHydroOnly",
    chi_radhydro::lua_utils::chiRadHydroMakePostShockConditionsHydroOnly);
  chi::console.RegisterFunction(
    "chiRadHydroTestMFEM",
    chi_radhydro::lua_utils::chiRadHydroTestMFEM);

  chi::RunBatch(argc,argv);

  chi::Finalize();

  chi::log.Log() << "chiRadHydro - Execution finished";

  return 0;
}