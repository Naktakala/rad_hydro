#include "rhsolverC.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

void chi_radhydro::SolverC_SNCN_MFEM::Initialize()
{
  SolverA_GDCN::Initialize();

  m_pwlc = SDM_PWLC::New(m_grid,
                         chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
                         chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  m_pwld = SDM_PWLD::New(m_grid,
                         chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
                         chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  const auto& grid = *m_grid;

  const auto ONE_DOF_PER_NODE =
    chi_math::UnknownManager({{chi_math::UnknownType::SCALAR,0}});

  //============================================= Prepare coefficient forms
  std::vector<MatDbl> list_of_Cc_star(grid.local_cells.size());
  std::vector<MatDbl> list_of_Cvol_star(grid.local_cells.size());

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c        = cell.local_id;
    const auto&    fem_data = m_pwlc->GetUnitIntegrals(cell);
    const size_t   Nd       = 3;
    const size_t   Nn       = fem_data.NumNodes();
    const size_t   Nf       = cell.faces.size();

    const auto& G = fem_data.GetIntV_gradshapeI();
    const auto& M = fem_data.GetIntV_shapeI_shapeJ();
    const auto& M_surf = fem_data.GetIntS_shapeI_shapeJ();
    const auto& V_c = fem_data.GetIntV_shapeI();

    //====================================== Build C-matrix
    MatDbl Cc(Nd * Nn, VecDbl(Nn+1, 0.0));
    //Left most column
    for (size_t i=0; i<Nn; ++i)
      for (size_t d=0; d<Nd; ++d)
        Cc[i*Nd + d][0] = - G[i][d];

    //Other columns
    for (size_t f=0; f<Nf; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& Mf = M_surf[f];
      const auto& nf = face.normal;

      for (size_t i=0; i<Nn; ++i)
        for (size_t j=0; j<Nn; ++j)
          for (size_t d=0; d<Nd; ++d)
            Cc[i*Nd + d][j+1] += nf[d] * Mf[i][j];
    }//for f

    //====================================== Build M-matrix
    MatDbl Mc(Nd*Nn, VecDbl(Nd*Nn, 0.0));
    for (size_t i=0; i<Nn; ++i)
      for (size_t j=0; j<Nn; ++j)
        for (size_t d=0; d<Nd; ++d)
          Mc[i*Nd + d][j*Nd + d] = M[i][j];

    //====================================== Build C_vol
    MatDbl Cvol(Nd * Nn, VecDbl(Nd, 0.0));
    for (size_t i=0; i<Nn; ++i)
      for (size_t d=0; d<Nd; ++d)
        Cvol[i*Nd+d][d] = V_c[i];

    //====================================== Compute C_c* and C_vol*
    const auto Mc_inv = chi_math::Inverse(Mc);
    auto& Cc_star = list_of_Cc_star[c];
    auto& Cvol_star = list_of_Cvol_star[c];
    Cc_star = chi_math::MatMul(Mc_inv, Cc);
    Cvol_star = chi_math::MatMul(Mc_inv, Cvol);
  }//for cell

  m_list_of_Cc_star = std::move(list_of_Cc_star);
  m_list_of_Cvol_star = std::move(list_of_Cvol_star);

  m_num_local_cfem_nodes = m_pwlc->GetNumLocalDOFs(ONE_DOF_PER_NODE);
  m_num_local_dfem_nodes = m_pwld->GetNumLocalDOFs(ONE_DOF_PER_NODE);
}