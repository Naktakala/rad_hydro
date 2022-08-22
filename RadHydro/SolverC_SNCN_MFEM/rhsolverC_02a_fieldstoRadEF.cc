#include "rhsolverC.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

void chi_radhydro::SolverC_SNCN_MFEM::
FieldsToRadEF(const std::vector<UVector>& U_n,
              VecDbl& rad_E_n,
              std::vector<Vec3>& rad_F_n,
              std::vector<Vec3>& rad_F0_n)
{
  const double c_spd = speed_of_light_cmpsh;
  const size_t num_local_nodes = m_grid->local_cells.size();

  auto& temp_rad_E = m_scalar_fields.at("radE");
  for (size_t i=0; i<num_local_nodes; ++i)
    rad_E_n[i] = temp_rad_E[i];

  //======================= Compute nodal averaged radE
  // These are in the block below
  // the cell averaged
  std::vector<double> nodal_vols(m_num_local_cfem_nodes,0.0);
  const int64_t Nc = static_cast<int64_t>(num_local_nodes);
  for (const auto& cell : m_grid->local_cells)
  {
    const uint64_t c = cell.local_id;
    const auto& fem_data = m_pwlc->GetUnitIntegrals(cell);
    size_t f = 0;
    for (const auto& face : cell.faces)
    {
      for (int fj=0; fj<face.vertex_ids.size(); ++fj)
      {
        const int j = fem_data.FaceDofMapping(f,fj);
        const int64_t jmap = m_pwlc->MapDOFLocal(cell, j);

        nodal_vols[jmap] += fem_data.IntV_shapeI(j);

        rad_E_n[jmap+Nc] += fem_data.IntV_shapeI(j) * rad_E_n[c];

      }//for fj
      ++f;
    }//for face
  }//for cell

  for (size_t n=0; n<m_num_local_cfem_nodes; ++n)
    rad_E_n[n + Nc] /= nodal_vols[n];

  //======================= Compute F^n
//  const size_t Nd = 3;
//  for (const auto& cell : m_grid->local_cells)
//  {
//    const uint64_t c = cell.local_id;
//    const int mat_id = cell.material_id;
//    const auto& fem_data = m_pwlc->GetUnitIntegrals(cell);
//    const double T = IdealGasTemperatureFromCellU(U_n[c],m_Cv);
//    const double kappa_a = ComputeKappaFromLua(T,mat_id, m_kappa_a_function);
//    const double kappa_s = ComputeKappaFromLua(T,mat_id, m_kappa_s_function);
//    const double kappa_t = kappa_a + kappa_s;
//
//    const double rho = U_n[c][RHO];
//    const Vec3   u   = VelocityFromCellU(U_n[c]);
//    const double sigma_t = rho*kappa_t;
//
//    const double D = c_spd / sigma_t;
//
//    const double VEFf = 1.0 / 3.0;
//
//    const auto& Cc_star = m_list_of_Cc_star[c];
//    const auto& Cvol_star = m_list_of_Cvol_star[c];
//
//    //Build colum map
//    const size_t Nn = fem_data.NumNodes();
//    std::vector<int64_t> Cc_star_col_map4knowns(Nn + 1, 0);
//    Cc_star_col_map4knowns[0] = static_cast<int64_t>(c);
//    for (size_t n=0; n<Nn; ++n)
//    {
//      const int64_t nmap = m_pwlc->MapDOFLocal(cell, n) + static_cast<int64_t>(Nc);
//      Cc_star_col_map4knowns[n + 1] = nmap;
//    }
//
//    const auto& known_col_map = Cc_star_col_map4knowns;
//    for (size_t i=0; i<Nn; ++i)
//    {
//      const int64_t imap = m_pwld->MapDOFLocal(cell, i);
//      auto& radF_c_i_n  = rad_F_n[imap];
//      auto& radF0_c_i_n = rad_F0_n[imap];
//
//      for (size_t d=0; d<Nd; ++d)
//        for (size_t n=0; n<(Nn + 1); ++n)
//          radF_c_i_n(d) += -D * Cc_star[i * Nd + d][n] *
//                           VEFf * rad_E_n[known_col_map[n]];
//
//      radF0_c_i_n = radF_c_i_n;
//
//      for (size_t d=0; d<Nd; ++d)
//        radF_c_i_n(d) +=
//          Cvol_star[i * Nd + d][d] *
//          (1.0 + VEFf) * rad_E_n[c] * u[d];
//
//    }//for i
//  }//for cell
}