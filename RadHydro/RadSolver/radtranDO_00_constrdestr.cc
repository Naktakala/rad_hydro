#include "radtranDO.h"

chi_radhydro::RadTranDO::RadTranDO(const std::string& name) :
  chi_physics::Solver(name,
                      {{"num_angles"      , int64_t(8)},
                       {"scattering_order", int64_t(0)},
                       {"num_groups"      , int64_t(0)}})
{}