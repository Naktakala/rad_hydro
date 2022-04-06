#include "compinfflow.h"

chi_hydro::CompInFFlow::
  CompInFFlow(const std::string& text_name) :
  chi_physics::Solver(text_name, {{"maximum_dt", double(1.0e-8)},
                                  {"CFL"       , double(0.9)   }})
{
  //nothing
}