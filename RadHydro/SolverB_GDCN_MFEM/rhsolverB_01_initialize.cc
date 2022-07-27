#include "rhsolverB.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

void chi_radhydro::SolverB_GDCN_MFEM::Initialize()
{
  SolverA_GDCN::Initialize();

  m_pwl = SDM_PWLC::New(m_grid,
                      chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
                      chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  const auto& grid = *m_grid;

  const auto ONE_DOF_PER_NODE =
    chi_math::UnknownManager({{chi_math::UnknownType::SCALAR,0}});

  //============================================= Prepare coefficient forms
  std::vector<MatDbl> list_of_Cc_star(grid.local_cells.size());

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c        = cell.local_id;
    const double   Dc       = 1.0;
    const auto&    fem_data = m_pwl->GetUnitIntegrals(cell);
    const size_t   Nd       = 3;
    const size_t   Nn       = fem_data.NumNodes();
    const size_t   Nf       = cell.faces.size();

    const auto& G = fem_data.GetIntV_gradshapeI();
    const auto& M = fem_data.GetIntV_shapeI_shapeJ();
    const auto& M_surf = fem_data.GetIntS_shapeI_shapeJ();

    //====================================== Build C-matrix
    MatDbl Cc(Nd * Nn, VecDbl(Nn+1, 0.0));
    //Left most column
    for (size_t i=0; i<Nn; ++i)
      for (size_t d=0; d<Nd; ++d)
        Cc[i*Nd + d][0] = -Dc * G[i][d];

    //Other columns
    for (size_t f=0; f<Nf; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& Mf = M_surf[f];
      const auto& nf = face.normal;

      for (size_t i=0; i<Nn; ++i)
        for (size_t j=0; j<Nn; ++j)
          for (size_t d=0; d<Nd; ++d)
            Cc[i*Nd + d][j+1] += Dc * nf[d] * Mf[i][j];
    }//for f

    //====================================== Build M-matrix
    MatDbl Mc(Nd*Nn, VecDbl(Nd*Nn, 0.0));
    for (size_t i=0; i<Nn; ++i)
      for (size_t j=0; j<Nn; ++j)
        for (size_t d=0; d<Nd; ++d)
          Mc[i*Nd + d][j*Nd + d] = M[i][j];

    //====================================== Compute C*
    const auto Mc_inv = chi_math::Inverse(Mc);
    auto& Cc_star = list_of_Cc_star[c];
    Cc_star = chi_math::MatMul(Mc_inv, Cc);
  }//for cell

  m_list_of_Cc_star = std::move(list_of_Cc_star);

  m_num_local_cfem_nodes = m_pwl->GetNumLocalDOFs(ONE_DOF_PER_NODE);
}