#include "lua_solver_base.h"

namespace chi_radhydro::solver_base_lua_utils
{

//###################################################################
/**Prescribes a bc.
\param solver_handle int    Handle to the relevant solver.
\param boundary_id   int    Boundary-id
\param boundary_type int    0 for Transmissive and 1 for fixed.
\param rho           double Density.
\param u             double Velocity x-component.
\param v             double Velocity y-component.
\param w             double Velocity z-component.
\param e             double Internal energy.
\param p             double Pressure.
\param radE          double Radiation energy.
*/
int chiRadHydroSetBCSetting(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  //============================== Get solver
  const int solver_handle = lua_tointeger(L, 1);

  auto& rad_hydro_solver = chi::GetStackItem<RadHydroSolver>(
    chi::solver_stack, solver_handle);

  const uint64_t boundary_id   = lua_tointeger(L, 2);
  const int      boundary_type = lua_tointeger(L, 3);

  if (boundary_type == 1 and num_args < 10)
    LuaPostArgAmountError(fname, 10, num_args);

  BCType bc_type_enum = BCType::TRANSMISSIVE;
  std::array<double, 7> bvalues = {0,0,0,0,0,0,0};

  if (boundary_type == 1)
  {
    bc_type_enum = BCType::FIXED;
    for (int k=0; k<7; ++k)
      bvalues[k] = lua_tonumber(L,k+4);
  }

  std::stringstream bvalues_str;
  for (int k=0; k<7; ++k)
    bvalues_str << bvalues[k] << " ";

  BCSetting setting = {bc_type_enum, bvalues};

  rad_hydro_solver.m_bc_settings[boundary_id] = setting;

  chi::log.Log() << "RadHydroSolver: BC set for id " << boundary_id
                 << " type " << boundary_type
                 << " values " << bvalues_str.str();

  return 0;
}

}//namespace chi_radhydro::solver_base_lua_utils