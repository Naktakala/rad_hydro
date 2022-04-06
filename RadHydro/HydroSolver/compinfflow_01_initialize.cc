#include "compinfflow.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics& chi_physics_handler;

#include "ChiMath/chi_math.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

//###################################################################
/**Initializes the solver.*/
void chi_hydro::CompInFFlow::Initialize()
{
  chi_log.Log() << "\nInitializing CompInFFlow solver";

  //======================================== Parse options
  delta_t_max = basic_options("maximum_dt").FloatValue();

  chi_log.Log() << "  CompInFFlow: delta_t_max set to "
                << std::setprecision(4) << std::scientific
                << delta_t_max;

  //======================================== Check grid
  if (regions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "  CompInFFlow: No regions added to solver.";
    exit(EXIT_FAILURE);
  }
  grid = regions.back()->GetGrid();

  if (grid == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << " CompInFFlow: No grid available from region.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check materials
  {
    size_t num_invalid_mat_cells = 0;
    for (const auto& cell : grid->local_cells)
    {
      if (cell.material_id < 0 or
          cell.material_id >= chi_physics_handler.material_stack.size())
        ++num_invalid_mat_cells;
    }

    if (num_invalid_mat_cells > 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << " CompInFFlow: A total of " << num_invalid_mat_cells
        << " have invalid material ids.";
      exit(EXIT_FAILURE);
    }
  }

  //======================================== Determine geometry type
  if (grid->local_cells[0].Type() == chi_mesh::CellType::SLAB)
    coordinate_system = CoordinateSystem::ONED_SLAB;
  else if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYGON)
    coordinate_system = CoordinateSystem::TWOD_CARTESIAN;
  else if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYHEDRON)
    coordinate_system = CoordinateSystem::THREED_CARTESIAN;

  //======================================== Initialize unknown manager
  {
    using namespace chi_math;
    for (int i=1; i<=5; ++i) uk_man_U.AddUnknown(UnknownType::SCALAR);

    uk_man_U.SetUnknownTextName(0, "rho");
    uk_man_U.SetUnknownTextName(1, "rho_u");
    uk_man_U.SetUnknownTextName(2, "rho_v");
    uk_man_U.SetUnknownTextName(3, "rho_w");
    uk_man_U.SetUnknownTextName(4, "E");
  }

  //======================================== Initialize spatial discretizations
  fv = SpatialDiscretization_FV::New(grid);

  //======================================== Initialize vectors
  num_nodes_local = fv->GetNumLocalDOFs(ChiMath::UNITARY_UNKNOWN_MANAGER);
  num_dofs_local = fv->GetNumLocalDOFs(uk_man_U);

  U_old.resize(num_dofs_local);
  U_new = U_old;

  //======================================== Initialize material properties
  gamma.assign(num_nodes_local, 1.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& material = chi_physics_handler.material_stack[cell.material_id];

    for (const auto& property : material->properties)
      if (property->property_name == "Gamma")
      {
        const double gamma_value = property->GetScalarValue();
        if (gamma_value > 1.0)
          gamma[cell.local_id] = property->GetScalarValue();
        else
        {
          chi_log.Log(LOG_ALLERROR)
            << " CompInFFlow: Material " << cell.material_id
            << " has an invalid gamma property. Value must be >1.0 but is "
            << gamma_value << ".";
          exit(EXIT_FAILURE);
        }
      }
  }

  //======================================== Initialize fields
  rho.assign(num_nodes_local, 1.0);
  u.assign(num_nodes_local, 0.0);
  v=w=u;
  p.assign(num_nodes_local, 1.0);
  T.assign(num_nodes_local, 1.0);
  e.assign(num_nodes_local, 1.0);
  Cv.assign(num_nodes_local, 1.0);

  //======================================== Apply field settings
//  bool rho_specified = false;
//  bool p_specified = false;
//  bool T_specified = false;
  bool e_specified = false;
  for (const auto& field_setting : logvol_field_settings)
  {
    const auto logvol      = std::get<0>(field_setting);
    const auto field_name  = std::get<1>(field_setting);
    const auto field_value = std::get<2>(field_setting);

    for (const auto& cell : grid->local_cells)
      if (logvol->Inside(cell.centroid))
      {
        if (field_name == "rho") {rho[cell.local_id] == field_value; }
        if (field_name == "u"  )  u  [cell.local_id] == field_value;
        if (field_name == "v"  )  v  [cell.local_id] == field_value;
        if (field_name == "w"  )  w  [cell.local_id] == field_value;
        if (field_name == "p"  ) {p  [cell.local_id] == field_value; }
        if (field_name == "T"  ) {T  [cell.local_id] == field_value; }
        if (field_name == "e"  ) {e  [cell.local_id] == field_value; e_specified = true;}
      }
  }//for field_setting

  if (not e_specified)
    for (uint64_t c=0; c<num_nodes_local; ++c)
      e[c] = p[c]/(rho[c]*(gamma[c]-1.0));

  for (uint64_t c=0; c<num_nodes_local; ++c)
    T[c] = e[c]/Cv[c];

  //======================================== Create field functions
  

  chi_log.Log() << "Done initializing CompInFFlow solver\n\n";
}