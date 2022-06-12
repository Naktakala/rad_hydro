#include "radtranGreyDiff.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"

void chi_radhydro::RadTranGreyDiffusion::Initialize()
{
  chi::log.Log() << "\nInitializing" + this->TextName() + "solver";

  //======================================== Parse options
  delta_t_max = basic_options("maximum_dt").FloatValue();

  chi::log.Log() << "  CompInFFlow: delta_t_max set to "
                 << std::setprecision(4) << std::scientific
                 << delta_t_max;

  C_cfl = basic_options("CFL").FloatValue();

  chi::log.Log() << "  CompInFFlow: C_cfl set to "
                 << std::setprecision(4) << std::scientific
                 << C_cfl;

  num_timesteps = basic_options("max_timesteps").IntegerValue();

  chi::log.Log() << "  CompInFFlow: num_timesteps set to "
                 << std::setprecision(4) << std::scientific
                 << num_timesteps;

  t_max = basic_options("max_time").FloatValue();

  chi::log.Log() << "  CompInFFlow: t_max set to "
                 << std::setprecision(4) << std::scientific
                 << t_max;

  //======================================== Check grid
  grid = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid == nullptr)
  {
    chi::log.LogAllError()
      << this->TextName() + ": No grid available from region.";
    chi::Exit(EXIT_FAILURE);
  }

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
  }

  //======================================== Initialize spatial discretizations
  fv = chi_math::SpatialDiscretization_FV::New(grid);

  //======================================== Initialize vectors
  auto unitary_uk_man = chi_math::UnknownManager::GetUnitaryUnknownManager();
  num_nodes_local = fv->GetNumLocalDOFs(unitary_uk_man);

  e_E_uk_man.AddUnknown(chi_math::UnknownType::SCALAR);
  e_E_uk_man.AddUnknown(chi_math::UnknownType::SCALAR);

  //======================================== Initialize material properties
  gamma.assign(num_nodes_local, 1.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& material = chi::GetStackItemPtr(chi::material_stack,cell.material_id);

    for (const auto& property : material->properties)
      if (property->property_name == "Gamma")
      {
        const double gamma_value = property->GetScalarValue();
        if (gamma_value > 1.0)
          gamma[cell.local_id] = property->GetScalarValue();
        else
        {
          chi::log.LogAllError()
            << " CompInFFlow: Material " << cell.material_id
            << " has an invalid gamma property. Value must be >1.0 but is "
            << gamma_value << ".";
          exit(EXIT_FAILURE);
        }
      }
  }

  //======================================== Initialize cell characteristic
  //                                         lenghts
  cell_char_length.assign(num_nodes_local, 0.0);
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

    cell_char_length[cell.local_id] = char_length;
  }

  //======================================== Initialize fields
  rho.assign(num_nodes_local, 1.0);
  u.assign(num_nodes_local, 0.0);
  v=w=u;
  p.assign(num_nodes_local, 1.0);
  temperature.assign(num_nodes_local, 1.0);
  e.assign(num_nodes_local, 1.0);
  rad_E = e;
  Cv.assign(num_nodes_local, 1.0);

  //======================================== Apply field settings
  bool e_specified = false;
  for (const auto& field_setting : logvol_field_settings)
  {
    const auto logvol      = std::get<0>(field_setting);
    const auto field_name  = std::get<1>(field_setting);
    const auto field_value = std::get<2>(field_setting);

    for (const auto& cell : grid->local_cells)
      if (logvol->Inside(cell.centroid))
      {
        if (field_name == "rho") {rho[cell.local_id] = field_value; }
        if (field_name == "u"  )  u  [cell.local_id] = field_value;
        if (field_name == "v"  )  v  [cell.local_id] = field_value;
        if (field_name == "w"  )  w  [cell.local_id] = field_value;
        if (field_name == "p"  ) {p  [cell.local_id] = field_value; }
        if (field_name == "T"  ) { temperature  [cell.local_id] = field_value; }
        if (field_name == "e"  ) {e  [cell.local_id] = field_value; e_specified = true;}
        if (field_name == "rad_E"  )  rad_E[cell.local_id] = field_value;
      }
  }//for field_setting

  if (not e_specified)
    for (uint64_t c=0; c<num_nodes_local; ++c)
      e[c] = p[c]/(rho[c]*(gamma[c]-1.0));

  for (uint64_t c=0; c<num_nodes_local; ++c)
    temperature[c] = e[c] / Cv[c];

  U_old.resize(num_nodes_local);
  FieldsToU(U_old, *grid, rho,u,v,w,e);

  //======================================== Create field functions
  if (field_functions.empty())
  {
    auto MakeFF = [this](const std::string& field_name,
                         std::vector<double>* field)
    {
      chi_math::UnknownManager uk_man;
      uk_man.AddUnknown(chi_math::UnknownType::SCALAR);
      uk_man.SetUnknownTextName(0, field_name);
      auto discretization = std::dynamic_pointer_cast<chi_math::SpatialDiscretization>(fv);
      auto ff = std::make_shared<chi_physics::FieldFunction>(
        TextName() + "-" + field_name,
        discretization,
        field,
        uk_man,
        0,0);

      chi::fieldfunc_stack.push_back(ff);
      field_functions.push_back(ff);
    };

    MakeFF("rho", &rho);
    MakeFF("u", &u);
    MakeFF("v", &v);
    MakeFF("w", &w);
    MakeFF("p", &p);
    MakeFF("e", &e);
    MakeFF("rad_E", &rad_E);
    MakeFF("T", &temperature);
    MakeFF("gamma", &gamma);
    MakeFF("Cv", &Cv);
  }//if field functions not created

  chi::log.Log() << "\nDone initializing" + this->TextName() + "solver";
}