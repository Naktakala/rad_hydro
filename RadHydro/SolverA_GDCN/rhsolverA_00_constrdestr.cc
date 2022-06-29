#include "rhsolverA.h"

chi_radhydro::SolverA_GDCN::SolverA_GDCN(const std::string &name) :
  RadHydroSolver(name)
{
  //======================================== Add basic option
  basic_options.AddOption("kappa_s_function",
                          std::string("MaterialKappaSFunction"));
  basic_options.AddOption("kappa_a_function",
                          std::string("MaterialKappaAFunction"));
  //======================================== Add scalar fields
  //Material properties
  scalar_fields["gamma"] = {};
  scalar_fields["Cv"   ] = {};

  //Unknowns
  scalar_fields["rho"  ] = {};
  scalar_fields["u"    ] = {};
  scalar_fields["v"    ] = {};
  scalar_fields["w"    ] = {};
  scalar_fields["p"    ] = {};
  scalar_fields["radE" ] = {};

  //Derived units
  scalar_fields["e"    ] = {};
  scalar_fields["temperature"] = {};
}