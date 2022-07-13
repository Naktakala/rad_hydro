#include "rhsolverA.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"

void chi_radhydro::SolverA_GDCN::Initialize()
{
  RadHydroSolver::Initialize();

  const size_t num_nodes_local = grid->local_cells.size();

  //======================================== Initialize material properties
  //gamma and Cv
  {
    size_t num_gammas_not_found = 0;
    size_t num_Cvs_not_found = 0;
    for (const auto& cell : grid->local_cells)
    {
      const auto& material = chi::material_stack[cell.material_id];

      bool gamma_found = false;
      bool Cv_found = false;
      for (const auto& property : material->properties)
      {
        if (property->property_name == "gamma")
        {
          const double value = property->GetScalarValue();
          if (value > 1.0)
            m_gamma = value;
          else
          {
            chi::log.LogAllError()
              << TextName() + ": Material " << cell.material_id
              << " has an invalid gamma property. Value must be >1.0 but is "
              << value << ".";
            exit(EXIT_FAILURE);
          }
          gamma_found = true;
        }//Gamma
        if (property->property_name == "Cv")
        {
          const double value = property->GetScalarValue();
          if (value > 0.01)
            m_Cv = value;
          else
          {
            chi::log.LogAllError()
              << TextName() + ": Material " << cell.material_id
              << " has an invalid Cv property. Value must be >0.01 but is "
              << value << ".";
            exit(EXIT_FAILURE);
          }
          Cv_found = true;
        }//Cv
      }//for property

      if (gamma_found and Cv_found) break;

      if (not gamma_found) ++num_gammas_not_found;
      if (not Cv_found) ++num_Cvs_not_found;
    }//for cell

    if (num_gammas_not_found > 0)
    {
      chi::log.LogAllError()
        << TextName() + ": Not all cells could be assigned gamma values. "
                        " Make sure all materials have the \"gamma\" scalar "
                        "property.";
      exit(EXIT_FAILURE);
    }
    if (num_Cvs_not_found > 0)
    {
      chi::log.LogAllError()
        << TextName() + ": Not all cells could be assigned Cv values. "
                        " Make sure all materials have the \"Cv\" scalar "
                        "property.";
      exit(EXIT_FAILURE);
    }
  }//gamma and Cv

  //======================================== Identify bndry cells
  for (const auto& cell : grid->local_cells)
    for (const auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        m_bndry_cells_local_ids.push_back(cell.local_id);
        break;
      }

  //======================================== Apply field settings
  bool rho_specified = false;
  bool e_specified = false;
  bool p_specified = false;
  bool T_specified = false;
  for (const auto& field_setting : logvol_field_settings)
  {
    const auto logvol      = std::get<0>(field_setting);
    const auto field_name  = std::get<1>(field_setting);
    const auto field_value = std::get<2>(field_setting);

    if (field_name == "rho"  )         rho_specified = true;
    if (field_name == "e"  )           e_specified = true;
    if (field_name == "temperature"  ) T_specified = true;
    if (field_name == "p")             p_specified = true;

    auto& field_data = scalar_fields.at(field_name);

    for (const auto& cell : grid->local_cells)
      if (logvol->Inside(cell.centroid))
        field_data[cell.local_id] = field_value;
  }//for field_setting

  //======================================== Compute deriveable settings
  if (not rho_specified) throw std::logic_error("rho not specified.");

  if (not p_specified and not(T_specified or e_specified))
    throw std::logic_error("If p is not specified then either T or e must be"
                           "specified.");
  if (not e_specified and not(T_specified or p_specified))
    throw std::logic_error("If e is not specified then either T or p must be"
                           "specified.");
  if (not T_specified and not(e_specified or p_specified))
    throw std::logic_error("If T is not specified then either e or p must be"
                           "specified.");

  //============================== set p if not specified
  if (not p_specified and T_specified)
  {
    auto& p = scalar_fields.at("p");

    const auto& rho = scalar_fields.at("rho");
    const auto& T = scalar_fields.at("temperature");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      p[c] = (m_gamma - 1)*rho[c]*m_Cv*T[c];
  }

  if (not p_specified and e_specified)
  {
    auto& p = scalar_fields.at("p");

    const auto& rho = scalar_fields.at("rho");
    const auto& e = scalar_fields.at("e");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      p[c] = (m_gamma - 1)*rho[c]*e[c];
  }

  //============================== set e if not specified
  if (not e_specified and T_specified)
  {
    auto& e = scalar_fields.at("e");

    const auto& T = scalar_fields.at("temperature");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      e[c] = T[c]*m_Cv;
  }

  if (not e_specified and p_specified)
  {
    auto& e = scalar_fields.at("e");

    const auto& p = scalar_fields.at("p");
    const auto& rho = scalar_fields.at("rho");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      e[c] = p[c]/(rho[c]*(m_gamma-1.0));
  }

  //============================== set p if not specified
  if (not T_specified and e_specified)
  {
    auto& temperature = scalar_fields.at("temperature");

    const auto& e  = scalar_fields.at("e");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      temperature[c] = e[c] / m_Cv;
  }

  if (not T_specified and p_specified)
  {
    auto& temperature = scalar_fields.at("temperature");

    const auto& p  = scalar_fields.at("p");
    const auto& rho = scalar_fields.at("rho");

    for (uint64_t c=0; c<num_nodes_local; ++c)
      temperature[c] = p[c]/(m_gamma-1)/rho[c]/m_Cv;
  }

  {
    const auto& temperature = scalar_fields.at("temperature");
    auto& radE = scalar_fields.at("radE");
    for (uint64_t c=0; c<num_nodes_local; ++c)
      radE[c] = a*pow(temperature[c],4.0);
//      radE[c] = 0.0;
  }



}