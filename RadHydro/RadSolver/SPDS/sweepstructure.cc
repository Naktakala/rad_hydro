#include "sweepstructure.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiGraph/chi_directed_graph.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Computes the sweep ordering.*/
chi_radhydro::SweepStructure::
  SweepStructure(const chi_mesh::MeshContinuum& grid,
                 const chi_mesh::Vector3& input_omega,
                 double parallel_tolerance/* = 1.0e-12*/) :
  m_master_direction(input_omega),
  m_parallel_tolerance(parallel_tolerance)
{
  /**Checks if a cell is a linear type.*/
  auto IsLinear = [](const chi_mesh::CellType& type)
  {
    if (type == chi_mesh::CellType::SLAB or
        type == chi_mesh::CellType::POLYGON or
        type == chi_mesh::CellType::POLYHEDRON)
      return true;
    return false;
  };

  //======================================== Initialize directed graph
  const size_t num_local_cells = grid.local_cells.size();

  chi_graph::DirectedGraph local_DG;

  for (uint64_t c=0; c<num_local_cells; ++c)
    local_DG.AddVertex();

  //======================================== Build local connectivity and
  //                                         local directed graph
  m_connectivity.resize(num_local_cells);

  const auto& omega = m_master_direction;
  const double& epsilon = m_parallel_tolerance;
  for (const auto& cell : grid.local_cells)
  {
    std::vector<Orientation> face_point_orientations;

    ChiInvalidArgument(not IsLinear(cell.Type()),
      "Only linear cell types are supported.")

    const size_t num_face_points = cell.faces.size();
    face_point_orientations.reserve(num_face_points);

    for (const auto& face : cell.faces)
    {
      const auto& n = face.normal;
      const double mu = n.Dot(omega);
      const bool has_local_neighbor = face.has_neighbor;

      uint64_t neighbor_local_id;
      if (has_local_neighbor)
        neighbor_local_id = face.GetNeighborLocalID(grid);

      Orientation             orientation = Orientation::PARALLEL;
      if (mu > epsilon)       orientation = Orientation::OUTGOING;
      else if (mu < -epsilon) orientation = Orientation::INCOMING;

      if (orientation == Orientation::OUTGOING and has_local_neighbor)
        local_DG.AddEdge(static_cast<int>(cell.local_id),
                         static_cast<int>(neighbor_local_id));
      else if (orientation == Orientation::INCOMING and has_local_neighbor)
        local_DG.AddEdge(static_cast<int>(neighbor_local_id),
                         static_cast<int>(cell.local_id));

      face_point_orientations.push_back(orientation);
    }//for face

    m_connectivity.push_back(std::move(face_point_orientations));
  }//for cell

  //======================================== Remove cyclic dependencies
  auto delayed_edges = local_DG.RemoveCyclicDependencies();

  //======================================== Generate topological sorting
  m_ordering = local_DG.GenerateTopologicalSort();
}