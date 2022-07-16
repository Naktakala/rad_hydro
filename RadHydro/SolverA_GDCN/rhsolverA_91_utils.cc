#include "rhsolverA.h"

void chi_radhydro::SolverA_GDCN::
  ProcessAndExportFields(const std::string &file_name,
                         const std::vector<UVector> &U,
                         const std::vector<double> &rad_E,
                         const size_t num_local_nodes,
                         const std::map<std::string,double>& scalar_properties)
{
  //======================================== Parse U and radE to field-functions
  UToFields(U,
            m_scalar_fields.at("rho"),
            m_scalar_fields.at("u"),
            m_scalar_fields.at("v"),
            m_scalar_fields.at("w"),
            m_scalar_fields.at("e"),
            m_scalar_fields.at("p"),
            m_gamma);

  m_scalar_fields.at("radE") = rad_E;

  //======================================== Compute material temperature
  auto& temperature = m_scalar_fields.at("temperature");

  const auto& e  = m_scalar_fields.at("e");

  for (uint64_t c=0; c < num_local_nodes; ++c)
    temperature[c] = e[c] / m_Cv;

  //======================================== Print outputs
  PrintRawOutput(file_name, scalar_properties);
}