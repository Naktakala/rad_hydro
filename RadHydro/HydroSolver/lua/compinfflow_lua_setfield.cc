#include "compinfflow_lua.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

namespace chi_hydro::compinfflow_lua_utils
{

//###################################################################
/**Passes a setting of a flow-field variable. The setting is done with a
 * logical volume. If a cell is within the logical volume the indicated field
 * value will be set.
\param solver_handle int         Handle to the relevant solver.
\param logical_volume_handle int Handle to the logical volume used to set the
                                 field value.
\param field_name string         String name of the field to set.
\param field_value double        Value for the field.

## _

### Field names
The following field names are supported:
 - "rho" The fluid density.
 - "u", "v", "w", The fluid x,y,z velocity respectively.
 - "p" The fluid pressure
 - "T" The fluid temperature

*/
int chiCompInFFlowSetFieldInitialValue(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError(fname, 4, num_args);

  for (int i=1; i<=4; ++i) LuaCheckNilValue(fname, L, 1);

  //============================== Get solver
  const int solver_handle = lua_tointeger(L, 1);

  auto& flow_solver = chi_hydro::compinfflow_lua_utils::
    GetSolverByHandle(solver_handle, fname);

  //============================== Get logical volume
  const int logvol_handle = lua_tointeger(L, 2);

  auto logvol = chi::GetStackItemPtr(chi::logicvolume_stack,
                                     logvol_handle, fname);

  //============================== Get name and value
  LuaCheckStringValue(fname, L, 3);
  const std::string field_name = lua_tostring(L, 3);
  const double      field_value = lua_tonumber(L, 4);

  //============================== Check names and values
  /**Lamda to check if string is in list.*/
  auto StringInList = [](
    const std::string& value, const std::vector<std::string>& list)
  {

    if (std::any_of(list.begin(), list.end(),
                    [&value](const std::string& list_value)
                    {return value == list_value;}
                    ))
      return true;
    return false;
  };

  if (not StringInList(field_name, {"rho","u","v","w","p","T","e", "Cv"}))
  {
    throw std::invalid_argument(fname + ": Invalid field-name " + field_name +
    " specified. Allowed options are rho, u, v, w, p, T, e and Cv.");
  }

  if (StringInList(field_name, {"rho","p","T","e", "Cv"}) and field_value <= 0.0)
    throw std::invalid_argument(fname + ": Invalid field-value " +
    std::to_string(field_value) + " specified for field " + field_name +
    ". Expected value must be > 0.0");

  //============================== Pass settings to solver
  flow_solver.logvol_field_settings.emplace_back(logvol, field_name, field_value);

  return 0;
}

}//namespace chi_hydro::compinfflow_lua_utils
