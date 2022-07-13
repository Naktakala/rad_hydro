#include "solver_base.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

namespace chi_radhydro
{

void RadHydroSolver::Initialize()
{
  chi::log.Log() << "Initializing RadHydroSolver entities.";

  //======================================== Check grid
  grid = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid == nullptr)
  {
    chi::log.LogAllError()
      << this->TextName() + ": No grid available. Make sure a grid is defined.";
    chi::Exit(EXIT_FAILURE);
  }

  //======================================== Initialize spatial discretizations
  fv = chi_math::SpatialDiscretization_FV::New(grid);

  //======================================== Check materials
  {
    size_t num_invalid_mat_cells = 0;
    for (const auto& cell : grid->local_cells)
    {
      if (cell.material_id < 0 or
          cell.material_id >= chi::material_stack.size())
        ++num_invalid_mat_cells;
    }

    if (num_invalid_mat_cells > 0)
    {
      chi::log.LogAllError()
        << this->TextName() + ": A total of " << num_invalid_mat_cells
        << " have invalid material ids.";
      chi::Exit(EXIT_FAILURE);
    }
  }//Check materials

  //======================================== Initialize scalar fields set
  //                                         during construction
  const size_t num_local_cells = grid->local_cells.size();
  for (auto& field_name_vec : scalar_fields)
    field_name_vec.second.assign(num_local_cells, 0.0);

  //======================================== Scalar field for cell_char_length
  chi::log.Log() << "Computing cell characteristic lengths.";
  auto& CCL = scalar_fields.at("cell_char_length");
  {
    CCL.assign(num_local_cells,0.0);
    for (const auto& cell : grid->local_cells)
    {
      double char_length = 0.0;
      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        for (size_t fp=0; fp < num_faces; ++fp) //fp = f_prime
        {
          if (f == fp) continue;
          const auto& fc  = cell.faces[f ].centroid;
          const auto& fpc = cell.faces[fp].centroid;

          const double L = (fpc - fc).Norm();
          char_length = std::max(L, char_length);
        }//fp
      }//f

      CCL[cell.local_id] = char_length;
    }
  }//cell_char_length

  //======================================== Create field functions
  chi::log.Log() << "Creating field functions.";
  if (field_functions.empty())
  {
    /**Lambda to make a FV field function*/
    auto MakeFF = [this](const std::string& field_name,
                         std::vector<double>* field)
    {
      using SD = chi_math::SpatialDiscretization;
      chi_math::UnknownManager uk_man;
      uk_man.AddUnknown(chi_math::UnknownType::SCALAR);
      uk_man.SetUnknownTextName(0, field_name);
      auto discretization = std::dynamic_pointer_cast<SD>(fv);
      auto ff = std::make_shared<chi_physics::FieldFunction>(
        TextName() + "-" + field_name,
        discretization,
        field,
        uk_man,
        0,0);

      chi::fieldfunc_stack.push_back(ff);
      field_functions.push_back(ff);
    };

    for (auto& field_name_vec : scalar_fields)
    {
      auto& field_name = field_name_vec.first;
      auto& field_data  =field_name_vec.second;
      MakeFF(field_name, &field_data);
    }
  }//if field functions not created

  //======================================== Initialize boundary conditions
  std::set<int> unique_boundary_ids;
  for (const auto& cell : grid->local_cells)
    for (const auto& face : cell.faces)
      if (not face.has_neighbor)
        unique_boundary_ids.insert(static_cast<int>(face.neighbor_id));

  for (int bid : unique_boundary_ids)
    if (bc_settings.count(bid) == 0)
      throw std::logic_error("Boundary with id " + std::to_string(bid) +
      " not specified.");

  chi::log.Log() << "Done initializing RadHydroSolver entities.";

}//Initialize

}//namespace chi_radhydro