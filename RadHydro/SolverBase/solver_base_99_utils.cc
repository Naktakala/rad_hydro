#include "solver_base.h"

#include "ChiTimer/chi_timer.h"

#include <fstream>

namespace chi_radhydro
{

//###################################################################
/**Prints all of the scalar fields to a text file.
 * \param file_name Name of the text file to write.
 * \param scalar_properties A list of scalar-values to print on the
 *         first line of the file. i.e. time=0.01 balance=4.2e-6 */
void RadHydroSolver::
  PrintRawOutput(const std::string &file_name,
                 const std::map<std::string,double>& scalar_properties)
{
  chi::log.Log() << TextName() << ": Writing raw output to "
                 << "\"" << file_name << "\"";
  const size_t num_local_nodes = m_grid->local_cells.size();

  //============================================= Open file
  std::ofstream file(file_name, std::ios_base::trunc);

  if (not file.is_open())
  {
    chi::log.Log0Error() << __PRETTY_FUNCTION__
                         << ": Could not open file \"" << file_name << "\"";
    return;
  }
  //============================================= Write scalar properties
  // At minimum always write time
  if (scalar_properties.count("time") == 0)
  {
    file << "time -1.0 ";
    for (const auto& name_value : scalar_properties)
      file << name_value.first << " " << name_value.second << " ";
    file << "\n";
  }
  else
  {
    file << "time " << scalar_properties.at("time") << " ";
    for (const auto& name_value : scalar_properties)
      if (name_value.first != "time")
        file << name_value.first << " " << name_value.second << " ";
    file << "\n";
  }


  //============================================= Write field names
  file << "x y z ";
  for (const auto& fn_field : m_scalar_fields)
    file << fn_field.first << " ";
  file << "\n";

  //============================================= Write the fields
  for (size_t k=0; k<num_local_nodes; ++k)
  {
    const auto& xc = m_grid->local_cells[k].centroid;
    file << xc.x << " " << xc.y << " " << xc.z << " ";
    for (const auto& fn_field : m_scalar_fields)
      file << fn_field.second[k] << " ";
    file << "\n";
  }//for k

  file.close();
}

//###################################################################
/**Processes a string containing floating point values and converts
 * it to a vector of doubles.*/
std::vector<double> RadHydroSolver::
  MakeOutputTimesFromStr(const std::string &output_times_str)
{
  if (output_times_str.empty()) return {};

  std::vector<double> output_times;

  std::istringstream iss(output_times_str);
  while (iss.tellg() >= 0)
  {
    double value; iss >> value;
    if (not (iss.rdstate() & std::istringstream::failbit))
      output_times.push_back(value);
  }

  return output_times;
}



}//namespace chi_radhydro