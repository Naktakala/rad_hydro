#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "ChiConsole/chi_console.h"

//###################################################################
/**This function calls a lua function with two arguments. The first
 * argument is a scalar, and the second argument is an integer.*/
double chi_radhydro::
  ComputeKappaFromLua(double T, int mat_id, const std::string& lua_fname)
{
  if (lua_fname.empty()) return 0.0;

  double ret_val = 0.0;

  auto& L = chi::console.consoleState;

  int result = lua_getglobal(L,lua_fname.c_str());
  if (result != LUA_TFUNCTION)
  {
    lua_pop(L,1);
    throw std::logic_error(std::string(__FUNCTION__) + ": Tasked with "
    "calling lua-function \"" + lua_fname + "\" but the global definition was not"
    "found.");
  }
  lua_pushnumber(L,T);
  lua_pushinteger(L,mat_id);

  //2 arguments, 1 result, 0=original error object
  if (lua_pcall(L,2,1,0) == 0)
  {
    ret_val = lua_tonumber(L,-1);
  }
  lua_pop(L,1);

  return ret_val;
}

//###################################################################
/**Computes kappa_s and a, for the given U-state for each cell.*/
void chi_radhydro::
  ComputeCellKappas(const chi_mesh::MeshContinuum &grid,
                    const double               Cv,
                    const std::vector<UVector> &U,
                    const std::string &kappa_s_function,
                    const std::string &kappa_a_function,
                    std::vector<double> &kappa_a,
                    std::vector<double> &kappa_t)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;
    const int mat_id = cell.material_id;

    const double T_c = IdealGasTemperatureFromCellU(U[c], Cv);

    const double kappa_s_c = ComputeKappaFromLua(T_c, mat_id, kappa_s_function);
    const double kappa_a_c = ComputeKappaFromLua(T_c, mat_id, kappa_a_function);

    kappa_a[c] = kappa_a_c;
    kappa_t[c] = kappa_s_c + kappa_a_c;
  }//for cell
}