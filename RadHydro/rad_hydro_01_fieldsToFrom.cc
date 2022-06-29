#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Passes the fields-data to a cell U vectors.*/
void chi_radhydro::FieldsToU(std::vector<UVector> &pU,
                             const chi_mesh::MeshContinuum& grid,
                             const std::vector<double>& rho,
                             const std::vector<double>& u,
                             const std::vector<double>& v,
                             const std::vector<double>& w,
                             const std::vector<double>& e)
{
  if (pU.size() != grid.local_cells.size())
    throw std::logic_error("chi_radhydro::FieldsToU used with incompatible"
                           " U-vector dimension.");

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;
    pU[c](0) = rho[c];
    pU[c](1) = rho[c]*u[c];
    pU[c](2) = rho[c]*v[c];
    pU[c](3) = rho[c]*w[c];
    pU[c](4) = rho[c]*(0.5*u[c]*u[c] + 0.5*v[c]*v[c] + 0.5*w[c]*w[c] + e[c]);
  }//for cell
}

//###################################################################
/**Passes the fields-data to a cell U vectors.*/
void chi_radhydro::FieldsToU(const std::vector<double>& rho,
                             const std::vector<double>& u,
                             const std::vector<double>& v,
                             const std::vector<double>& w,
                             const std::vector<double>& e,
                             std::vector<UVector> &pU)
{
  const size_t num_nodes = pU.size();

  if (num_nodes != rho.size() or
      num_nodes != u.size() or
      num_nodes != v.size() or
      num_nodes != w.size() or
      num_nodes != e.size() )
    throw std::logic_error("chi_radhydro::FieldsToU used with incompatible"
                           " U-vector dimension.");

  for (size_t c=0; c<num_nodes; ++c)
  {
    pU[c](0) = rho[c];
    pU[c](1) = rho[c]*u[c];
    pU[c](2) = rho[c]*v[c];
    pU[c](3) = rho[c]*w[c];
    pU[c](4) = rho[c]*(0.5*u[c]*u[c] + 0.5*v[c]*v[c] + 0.5*w[c]*w[c] + e[c]);
  }//for cell
}

//###################################################################
/**Passes cell U vectors to fields data.*/
void chi_radhydro::
UToFields(const std::vector<UVector> &pU,
          const chi_mesh::MeshContinuum& grid,
          std::vector<double>& rho,
          std::vector<double>& u,
          std::vector<double>& v,
          std::vector<double>& w,
          std::vector<double>& e,
          std::vector<double>& p,
          const std::vector<double>& gamma)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;

    rho[c] = pU[c][0];
    u[c]   = pU[c][1]/rho[c];
    v[c]   = pU[c][2]/rho[c];
    w[c]   = pU[c][3]/rho[c];

    const double E = pU[c][4];
    e[c]   = E/rho[c] - 0.5*u[c]*u[c] - 0.5*v[c]*v[c] - 0.5*w[c]*w[c];
    p[c]   = (gamma[c]-1.0)*rho[c]*e[c];
  }//for cell
}

//###################################################################
/**Passes cell U vectors to fields data.*/
void chi_radhydro::
UToFields(const std::vector<UVector> &pU,
          std::vector<double>& rho,
          std::vector<double>& u,
          std::vector<double>& v,
          std::vector<double>& w,
          std::vector<double>& e,
          std::vector<double>& p,
          const std::vector<double>& gamma)
{
  const size_t num_nodes = pU.size();

  if (num_nodes != rho.size() or
      num_nodes != u.size() or
      num_nodes != v.size() or
      num_nodes != w.size() or
      num_nodes != e.size() or
      num_nodes != p.size() or
      num_nodes != gamma.size() )
    throw std::logic_error("chi_radhydro::FieldsToU used with incompatible"
                           " U-vector dimension.");

  for (size_t c=0; c<num_nodes; ++c)
  {
    rho[c] = pU[c][0];
    u[c]   = pU[c][1]/rho[c];
    v[c]   = pU[c][2]/rho[c];
    w[c]   = pU[c][3]/rho[c];

    const double E = pU[c][4];
    e[c]   = E/rho[c] - 0.5*u[c]*u[c] - 0.5*v[c]*v[c] - 0.5*w[c]*w[c];
    p[c]   = (gamma[c]-1.0)*rho[c]*e[c];
  }//for cell
}