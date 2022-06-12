#include "radtran_lua.h"

#include "../radtranDO.h"
#include "../radtranGreyDiff.h"
#include "chi_runtime.h"

//###################################################################
/**Creates a Discrete Ordinates Radiative Transfer solver.
\param solver_name string (Optional). Text name of the solver. [Default=RadTranDO]
\return solver_handle int
*/
int chi_radhydro::lua_utils::chiCreateRadTranDOSolver(lua_State *L)
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

//###################################################################
/**Creates a Grey Diffusion Radiative Transfer solver.
\param solver_name string (Optional). Text name of the solver. [Default=RadTranGreyDiff]
\return solver_handle int
*/
int chi_radhydro::lua_utils::chiCreateRadTranGreyDiffSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string solver_name = "RadTranGreyDiff";
  if (num_args >= 1)
  {
    LuaCheckNilValue(fname, L, 1);
    LuaCheckStringValue(fname, L, 1);

    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<RadTranGreyDiffusion>(solver_name);

  chi::solver_stack.push_back(new_solver);
  const size_t handle = chi::solver_stack.size()-1;
  lua_pushinteger(L, static_cast<lua_Integer>(handle));

  return 1;
}

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
int chi_radhydro::lua_utils::chiRadTranGreyDiffSetFieldInitialValue(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError(fname, 4, num_args);

  for (int i=1; i<=4; ++i) LuaCheckNilValue(fname, L, 1);

  //============================== Get solver
  const int solver_handle = lua_tointeger(L, 1);

  auto& flow_solver = chi::GetStackItem<chi_radhydro::RadTranGreyDiffusion>(
    chi::solver_stack, solver_handle);

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

void chi_radhydro::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lua_register(L, "chiCreateRadTranDOSolver", chiCreateRadTranDOSolver);
  lua_register(L, "chiCreateRadTranGreyDiffSolver", chiCreateRadTranGreyDiffSolver);
  lua_register(L, "chiRadTranGreyDiffSetFieldInitialValue", chiRadTranGreyDiffSetFieldInitialValue);
}