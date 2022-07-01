#include "solver_base.h"

#include <fstream>

namespace chi_radhydro
{

void RadHydroSolver::PrintRawOutput(const std::string &file_name)
{
  chi::log.Log() << TextName() << ": Writing raw output to "
                 << "\"" << file_name << "\"";
  const size_t num_local_nodes = grid->local_cells.size();

  //============================================= Open file
  std::ofstream file(file_name, std::ios_base::trunc);

  if (not file.is_open())
  {
    chi::log.Log0Error() << __PRETTY_FUNCTION__
                         << ": Could not open file \"" << file_name << "\"";
    return;
  }

  //============================================= Write field names
  file << "x y z ";
  for (const auto& fn_field : scalar_fields)
    file << fn_field.first << " ";
  file << "\n";

  //============================================= Write the fields
  for (size_t k=0; k<num_local_nodes; ++k)
  {
    const auto& xc = grid->local_cells[k].centroid;
    file << xc.x << " " << xc.y << " " << xc.z << " ";
    for (const auto& fn_field : scalar_fields)
      file << fn_field.second[k] << " ";
    file << "\n";
  }//for k

  file.close();
}

}//namespace chi_radhydro