#ifndef RADHYDRO_HYDRO_SOLVER_LUA_H
#define RADHYDRO_HYDRO_SOLVER_LUA_H

#include "chi_lua.h"

#include "../compinfflow.h"

namespace chi_hydro::compinfflow_lua_utils
{
  //###################################################################
  /** Obtains a pointer to a chi_hydro::CompInFFlow object or an object
   * derived from it.
   *
   * \param handle int Index in the chi_physics_handler where the solver object
   *                   should be located.
   * \param calling_function_name string The string used to print error messages,
   *                              should uniquely identify the calling function.
   *
   */
  chi_hydro::CompInFFlow*
  GetSolverByHandle(int handle, const std::string& calling_function_name);

  int chiCreateCompInFFlowSolver(lua_State* L);
  int chiCompInFFlowSetFieldInitialValue(lua_State* L);

  void RegisterLuaEntities(lua_State* L);
}//namespace chi_radhydro

#endif //RADHYDRO_HYDRO_SOLVER_LUA_H
