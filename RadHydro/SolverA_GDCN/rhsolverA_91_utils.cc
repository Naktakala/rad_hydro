#include "rhsolverA.h"

void chi_radhydro::SolverA_GDCN::
  ProcessAndExportFields(const std::string &file_name,
                         const std::vector<UVector> &U,
                         const std::vector<double> &rad_E,
                         const size_t num_local_nodes)
{
  //======================================== Parse U and radE to field-functions
  UToFields(U,
            scalar_fields.at("rho"),
            scalar_fields.at("u"),
            scalar_fields.at("v"),
            scalar_fields.at("w"),
            scalar_fields.at("e"),
            scalar_fields.at("p"),
            m_gamma);

  scalar_fields.at("radE") = rad_E;

  //======================================== Compute material temperature
  auto& temperature = scalar_fields.at("temperature");

  const auto& e  = scalar_fields.at("e");

  for (uint64_t c=0; c < num_local_nodes; ++c)
    temperature[c] = e[c] / m_Cv;

  //======================================== Print outputs
  PrintRawOutput(file_name);
}