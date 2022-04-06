#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "RadHydro/HydroSolver/lua/compinfflow_lua.h"

int main(int argc, char* argv[])
{
  ChiLog&     log     = ChiLog::GetInstance();

  log.Log(LOG_0) << "chiRadHydro - Execution started";

  ChiTech::Initialize(argc,argv);

  auto& lua_console = ChiConsole::GetInstance();

  chi_hydro::compinfflow_lua_utils::RegisterLuaEntities(lua_console.consoleState);

  ChiTech::RunBatch(argc,argv);

  ChiTech::Finalize();

  log.Log(LOG_0) << "chiRadHydro - Execution finished";

  return 0;
}