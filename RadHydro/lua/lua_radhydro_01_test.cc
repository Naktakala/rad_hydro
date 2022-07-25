#include "ChiLua/chi_lua.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_radhydro::lua_utils
{

int chiRadHydroTestMFEM(lua_State* L)
{
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;

  //============================================= Define lambda's
  auto PrintSystem = [](const MatDbl& mat, const VecDbl& rhs={},
                        const double div_fac_mat=1.0,
                        const double div_fac_rhs=1.0)
  {
    int j=0;
    for (const auto& row : mat)
    {
      for (double val : row)
        if (std::fabs(val) > 1.0e-12)
          printf("%5.1f ",val/div_fac_mat);
        else
          printf("00000 ");
      if (not rhs.empty())
        printf("     %5.1f ",rhs[j]/div_fac_rhs);
      std::cout << "\n";
      ++j;
    }
  };

  //============================================= Pickup stuff
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();

  const auto& grid = *grid_ptr;

  auto fv  = chi_math::SpatialDiscretization_FV::New(grid_ptr);
  auto pwl = chi_math::SpatialDiscretization_PWLC::New(grid_ptr,
                               chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
                               chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  const auto ONE_DOF_PER_NODE =
    chi_math::UnknownManager({{chi_math::UnknownType::SCALAR,0}});

  //============================================= Prepare coefficient forms
  chi::log.Log() << "Looping over cells.";
  std::vector<MatDbl> list_of_Cc(grid.local_cells.size());
  std::vector<MatDbl> list_of_Cc_star(grid.local_cells.size());

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c        = cell.local_id;
    const double   Dc       = -1.0 / 3.0;
    const auto&    fem_data = pwl->GetUnitIntegrals(cell);
    const size_t   Nd       = 3;
    const size_t   Nn       = fem_data.NumNodes();
    const size_t   Nf       = cell.faces.size();

    const auto& G = fem_data.GetIntV_gradshapeI();
    const auto& M = fem_data.GetIntV_shapeI_shapeJ();
    const auto& M_surf = fem_data.GetIntS_shapeI_shapeJ();

    //====================================== Compute volume
    double Vc = 0.0;
    for (size_t i=0; i<Nn; ++i) Vc += fem_data.IntV_shapeI(i);

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

//      if (face.has_neighbor)
      for (size_t i=0; i<Nn; ++i)
        for (size_t j=0; j<Nn; ++j)
          for (size_t d=0; d<Nd; ++d)
            Cc[i*Nd + d][j+1] += Dc * nf[d] * Mf[i][j];

    }//for f

    list_of_Cc[c] = Cc;

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

  //============================================= Develop cell node mapping
  // Suppose we have Nc number of cells. There would Nc number of unknowns for
  // RadE. The number of interstitial nodes would be hard to determine in any
  // coordinate system other than 1D. In 1D the number of interstitial nodes
  // would be Nc+1.
  const size_t num_radEc_nodes = fv->GetNumLocalDOFs(ONE_DOF_PER_NODE);
  const size_t num_radEj_nodes = pwl->GetNumLocalDOFs(ONE_DOF_PER_NODE);
  const size_t num_total_nodes = num_radEc_nodes + num_radEj_nodes;

  chi::log.Log() << "Number of internal nodes    : " << num_radEc_nodes;
  chi::log.Log() << "Number of interstitial nodes: " << num_radEj_nodes;

  MatDbl A(num_total_nodes, VecDbl(num_total_nodes, 0.0));
  VecDbl b(num_total_nodes, 0.0);

  //============================================= Build flags for boundary nodes
  std::vector<bool> node_bndry_flag(num_total_nodes, false);
  for (const auto& cell : grid.local_cells)
  {
    const auto& fem_data = pwl->GetUnitIntegrals(cell);
    const size_t num_faces = cell.faces.size();

    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];

      if (face.has_neighbor) continue;

      const size_t num_face_nodes = face.vertex_ids.size();
      for (size_t fj=0; fj<num_face_nodes; ++fj)
      {
        const size_t j = fem_data.FaceDofMapping(f, fj);
        const int64_t jmap = pwl->MapDOFLocal(cell, j) +
                             static_cast<int64_t>(num_radEc_nodes);
        node_bndry_flag[jmap] = true;
      }//for fj
    }//for face
  }

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c        = cell.local_id;
    const auto&    fem_data = pwl->GetUnitIntegrals(cell);
    const size_t   Nd       = 3;
    const size_t   Nn       = fem_data.NumNodes();
    const size_t   Nf       = cell.faces.size();

    const auto& S_surf = fem_data.GetIntS_shapeI();

    //====================================== Nodal map the columns of Cc_stars
    std::vector<int64_t> PCmap(Nn + 1, 0);
    PCmap[0] = fv->MapDOFLocal(cell, 0);
    for (size_t n=0; n<Nn; ++n)
      PCmap[n + 1] = pwl->MapDOFLocal(cell, n) +
                     static_cast<int64_t>(num_radEc_nodes);

    for (size_t j=0; j<Nn; ++j)
      b[c] += fem_data.IntV_shapeI(j);

    //====================================== Loop over faces
    const auto& Cc      = list_of_Cc[c];
    const auto& Cc_star = list_of_Cc_star[c];
    for (size_t f=0; f<Nf; ++f)
    {
      const auto&  face = cell.faces[f];
      const auto&  nf   = face.normal;
      const auto&  Sf   = S_surf[f];
      const size_t Nfn  = face.vertex_ids.size();

      //=============================== Assemble primary equation portion
      for (size_t fj=0; fj<Nfn; ++fj)
      {
        const size_t  j    = fem_data.FaceDofMapping(f, fj);

        for (size_t d = 0; d < Nd; ++d)
          for (size_t n = 0; n < (Nn + 1); ++n)
            if (node_bndry_flag[PCmap[n]])
              b[c] -= nf[d] * Sf[j] * Cc_star[j * Nd + d][n]*0.0;
            else
              A[c][PCmap[n]] += nf[d] * Sf[j] * Cc_star[j * Nd + d][n];
      }

      //=============================== Assembly auxiliary portions
      if (face.has_neighbor)
        for (size_t fj=0; fj<Nfn; ++fj)
        {
          const size_t  j    = fem_data.FaceDofMapping(f, fj);
          const int64_t jmap = pwl->MapDOFLocal(cell, j) +
                               static_cast<int64_t>(num_radEc_nodes);

          for (size_t d=0; d<Nd; ++d)
            for (size_t n=0; n<(Nn+1); ++n)
              if (node_bndry_flag[PCmap[n]])
                b[jmap] += nf[d] * Sf[j] * Cc_star[j * Nd + d][n]*0.0;
              else
                A[jmap][PCmap[n]] -= nf[d] * Sf[j] * Cc_star[j * Nd + d][n];
        }//for fj
      else
        for (size_t fj=0; fj<Nfn; ++fj)
        {
          const size_t  j    = fem_data.FaceDofMapping(f, fj);
          const int64_t jmap = pwl->MapDOFLocal(cell, j) +
                               static_cast<int64_t>(num_radEc_nodes);

          A[jmap][jmap] += 1.0;
          b[jmap] = 0.0;
        }//for fj
    }//for f

  }//for cell

  const double D = 1.0/3.0;
  const double h = 0.5/10;
  PrintSystem(A,b, 2*D/h, h);


  //Check symmetry
  for (size_t i=0; i<num_total_nodes; ++i)
    for (size_t j=0; j<num_total_nodes; ++j)
      if (std::fabs(A[i][j]-A[j][i]) > 1.0e-12)
      {
        chi::log.Log() << "Symmetry check failed.";
        break;
      }

  auto x = b;
  chi_math::GaussElimination(A,x, num_total_nodes);

//  x = chi_math::MatMul(A,x);

//  chi::log.Log() << "\n";
//  for (auto val : x)
//    printf("%.4e\n", val);

  for (size_t i=0; i<num_radEc_nodes; ++i)
  {
    printf("%.4e, %.4e\n", grid.cells[i].faces[0].centroid.z, x[i+num_radEc_nodes]);
    printf("%.4e, %.4e\n", grid.cells[i].centroid.z, x[i]);
  }
  const size_t N = num_radEc_nodes-1;
  printf("%.4e, %.4e\n", grid.cells[N].faces[1].centroid.z, x.back());

  return 0;
}

}//namespace chi_radhydro