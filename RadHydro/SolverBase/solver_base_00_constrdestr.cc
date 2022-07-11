#include "solver_base.h"

#include "ChiMath/UnknownManager/unknown_manager.h"

namespace chi_radhydro
{

//###################################################################
/**Constructor.*/
RadHydroSolver::RadHydroSolver(const std::string &name) :
  chi_physics::Solver(name),
  ONE_DOF_PER_NODE({{chi_math::UnknownType::SCALAR,1}})
{
  basic_options.AddOption("maximum_dt"    , double(1.0e-8)           );
  basic_options.AddOption("CFL"           , double(0.9)              );
  basic_options.AddOption("max_timesteps" , int64_t(1000)            );
  basic_options.AddOption("max_time"      , double(-1.0)             );
  basic_options.AddOption("export_times"  , std::string("")          );
  basic_options.AddOption("output_prefix" , std::string("ZRawOutput"));

  scalar_fields["cell_char_length"] = {};
}

}//namespace chi_radhydro