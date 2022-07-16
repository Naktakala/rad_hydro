#include "lua_solver_base.h"

namespace chi_radhydro::solver_base_lua_utils
{

//###################################################################
/**Passes a setting of a scalar variable. The setting is done with a
 * logical volume. If a cell is within the logical volume the indicated field
 * value will be set.
\param solver_handle int         Handle to the relevant solver.
\param logical_volume_handle int Handle to the logical volume used to set the
                                 field value.
\param field_name string         String name of the field to set.
\param field_value double        Value for the field.

*/
int chiRadHydroSolverSetScalarFieldWithLV(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError(fname, 4, num_args);

  for (int i=1; i<=4; ++i) LuaCheckNilValue(fname, L, 1);

  //============================== Get solver
  const int solver_handle = lua_tointeger(L, 1);

  auto& rad_hydro_solver = chi::GetStackItem<RadHydroSolver>(
    chi::solver_stack, solver_handle);

  //============================== Get logical volume
  const int logvol_handle = lua_tointeger(L, 2);

  auto logvol = chi::GetStackItemPtr(chi::logicvolume_stack,
                                     logvol_handle, fname);

  //============================== Get name and value
  LuaCheckStringValue(fname, L, 3);
  const std::string field_name = lua_tostring(L, 3);
  const double      field_value = lua_tonumber(L, 4);

  //============================== Pass settings to solver
  rad_hydro_solver.m_logvol_field_settings.emplace_back(
    logvol, field_name, field_value);

  return 0;
}

}//namespace chi_radhydro::solver_base_lua_utils