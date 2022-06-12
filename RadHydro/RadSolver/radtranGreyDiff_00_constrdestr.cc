#include "radtranGreyDiff.h"

chi_radhydro::RadTranGreyDiffusion::
  RadTranGreyDiffusion(const std::string &name) :
  chi_physics::Solver(name, {{"maximum_dt"   , double(1.0e-8)},
                             {"CFL"          , double(0.9)   },
                             {"max_timesteps", int64_t(1000) },
                             {"max_time"     , double(-1.0)}})
{}