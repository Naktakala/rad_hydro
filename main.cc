#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "RadHydro/HydroSolver/lua/compinfflow_lua.h"
#include "RadHydro/RadSolver/lua/radtran_lua.h"

int main(int argc, char* argv[])
{
  chi::log.Log() << "chiRadHydro - Execution started";

  chi::Initialize(argc,argv);

  chi_hydro::compinfflow_lua_utils::RegisterLuaEntities(chi::console.consoleState);
  chi_radtran::lua_utils::RegisterLuaEntities(chi::console.consoleState);

  chi::RunBatch(argc,argv);

  chi::Finalize();

  chi::log.Log() << "chiRadHydro - Execution finished";

  return 0;
}