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
  //Unknowns
  m_scalar_fields["rho"  ] = {};
  m_scalar_fields["u"    ] = {};
  m_scalar_fields["v"    ] = {};
  m_scalar_fields["w"    ] = {};
  m_scalar_fields["p"    ] = {};
  m_scalar_fields["radE" ] = {};

  //Derived units
  m_scalar_fields["e"    ] = {};
  m_scalar_fields["temperature"] = {};
}