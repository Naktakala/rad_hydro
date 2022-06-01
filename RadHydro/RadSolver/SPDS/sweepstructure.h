#ifndef RADHYDRO_SWEEPSTRUCTURE_H
#define RADHYDRO_SWEEPSTRUCTURE_H

#include <vector>

namespace chi_mesh
{
  class MeshContinuum;
}

#include "ChiMesh/chi_mesh.h"

namespace chi_radtran
{
  /**LinearElementSPDS*/
  class SweepStructure
  {
  public:
    enum class Orientation
    {
      INCOMING = 0,
      OUTGOING = 1,
      PARALLEL = 2
    };
  private:
    const chi_mesh::Vector3 m_master_direction;
    const double            m_parallel_tolerance;
    std::vector<std::vector<Orientation>> m_connectivity;
    std::vector<size_t> m_ordering;

  public:
    SweepStructure(const chi_mesh::MeshContinuum& grid,
                   const chi_mesh::Vector3& omega,
                   double parallel_tolerance = 1.0e-12);

  public:
    const chi_mesh::Vector3& MasterOmega() const {return m_master_direction;}
    double ParallelTolerance() const {return m_parallel_tolerance;}
    const std::vector<std::vector<Orientation>>& CellFacePointOrientations()
    {return m_connectivity;}
    const std::vector<size_t>& GridOrdering() const {return m_ordering;}
  };
}//namespace chi_radtran

#endif //RADHYDRO_SWEEPSTRUCTURE_H
