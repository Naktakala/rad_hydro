#include "rhsolverB.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "ChiMath/chi_math_banded_solvers.h"

void chi_radhydro::SolverB_GDCN_MFEM::
  AssySolveMFEMEnergyExchange(
  SimRefs&                        sim_refs,
  SDM_PWLC&                       pwl,
  std::vector<MatDbl>&            list_of_Cc_star,
  const std::vector<double>&      kappa_a_n,
  const std::vector<double>&      kappa_t_n,
  const std::vector<double>&      kappa_a_nph,
  const std::vector<double>&      kappa_t_nph,
  double                          tau,
  double                          theta1,
  double                          theta2,
  const std::vector<UVector>&     U_n,
  const std::vector<UVector>&     U_nph,
  const std::vector<UVector>&     U_nphstar,
  const std::vector<GradUTensor>& grad_U_nph,
  const std::vector<double>&      rad_E_n,
  const std::vector<double>&      rad_E_nph,
  const std::vector<double>&      rad_E_nphstar,
  const std::vector<Vec3>&        grad_rad_E_nph,
  std::vector<UVector>&           U_np1,
  std::vector<double>&            rad_E_np1)
{
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
          printf("%8.1e ",val/div_fac_mat);
        else
          printf("00000000 ");
      if (not rhs.empty())
        printf("     %5.1f ",rhs[j]/div_fac_rhs);
      std::cout << "\n";
      ++j;
    }
  };

  const auto&  grid = sim_refs.grid;
        auto&  fv   = sim_refs.fv;
  const double Cv   = sim_refs.Cv;

  const size_t Nc    = grid.local_cells.size(); //Number of cells
  const size_t Ns    = rad_E_n.size();
  const size_t Ncfem = Ns - Nc;

  MatDbl A;
  VecDbl b;
  A.assign(Ns, VecDbl(Ns, 0.0));
  b.assign(Ns, 0.0);

  sim_refs.e_recon_coeff_k5.assign(Ns, 0.0);
  sim_refs.e_recon_coeff_k6.assign(Ns, 0.0);

  const size_t Nd = 3; //Number of dimensions

  for (const auto& cell_c : grid.local_cells)
  {
    //=========================================== Get cell geometry info
    const uint64_t c         = cell_c.local_id;
    const size_t   Nf        = cell_c.faces.size();

    const auto&    fv_view_c   = fv.MapFeView(c);
    const auto&    face_areas  = fv_view_c->face_area;
    const Vec3&    x_c         = cell_c.centroid;
    const double   V_c         = fv_view_c->volume;

    //=========================================== Get fem stuff
    const auto&    fem_data = pwl.GetUnitIntegrals(cell_c);
    const size_t   Nn       = fem_data.NumNodes();

    const auto&    S_surf   = fem_data.GetIntS_shapeI();

    //=========================================== Get cell physics-info
    const double   rho_c_n     = U_n[c][RHO];
    const double   rho_c_np1   = U_np1[c][RHO];
    const double   T_c_n       = IdealGasTemperatureFromCellU(U_n[c], Cv);
    const double   T_c_nphstar = IdealGasTemperatureFromCellU(U_nphstar[c], Cv);

    const double   E_c_nphstar = U_nphstar[c][MAT_E];

    const double   e_c_nphstar     = InternalEnergyFromCellU(U_nphstar[c]);
    const double   rad_E_c_n       = rad_E_n[c];
    const double   rad_E_c_nphstar = rad_E_nphstar[c];

    const Vec3     u_c_np1         = VelocityFromCellU(U_np1[c]);

    //=========================================== Compute sigmas
    const double sigma_t_c_n   = rho_c_n * kappa_t_n[c];
    const double sigma_a_c_n   = rho_c_n * kappa_a_n[c];
    const double sigma_t_c_np1 = rho_c_np1 * kappa_t_nph[c];
    const double sigma_a_c_np1 = rho_c_np1 * kappa_a_nph[c];

    //=========================================== Catch infinite diffusion coeffs
    if (std::fabs(sigma_t_c_np1) < 1.0e-6)
    {
      A[c][c] = 1.0; b[c] = rad_E_c_nphstar;
      sim_refs.e_recon_coeff_k5[c] = 0.0;
      sim_refs.e_recon_coeff_k6[c] = e_c_nphstar;
      continue;
    }

    //=========================================== Compute explicit terms
    double third_grad_rad_E_dot_u_nph =
      Make3rdGradRadEDot_u(sim_refs,
                           cell_c, V_c, face_areas,
                           rad_E_nph, grad_rad_E_nph,
                           U_nph, grad_U_nph);

    // sigma_a c (aT^4 - radE)
    double Sea_n = MakeEmAbsSource(sigma_a_c_n, T_c_n, rad_E_c_n);

    //=========================================== Compute constants
    const double k1 = theta1 * sigma_a_c_np1 * speed_of_light_cmpsh;
    const double k2 = 4 * pow(T_c_nphstar, 3) / Cv;

    const double k3 = - k1*a*pow(T_c_nphstar, 4)
                      + k1*a*k2*e_c_nphstar
                      - theta2 * Sea_n
                      - third_grad_rad_E_dot_u_nph;
    const double k4 = -k1*a*k2;

    const double k5 = k1/(tau * rho_c_np1 - k4);
    const double k6 =
      (k3 - tau * 0.5 * rho_c_np1 * u_c_np1.NormSquare() + tau * E_c_nphstar) /
      (tau * rho_c_np1 - k4);

    sim_refs.e_recon_coeff_k5[c] = k5;
    sim_refs.e_recon_coeff_k6[c] = k6;

    //====================================== Nodal map the columns of Cc_stars
    //TODO: Begin : Find generalization of mapping
    std::vector<int64_t> Cc_star_col_map(Nn + 1, 0);

    Cc_star_col_map[0] = static_cast<int64_t>(2*c+1);
    for (size_t n=0; n<Nn; ++n)
      if (n==0)
        Cc_star_col_map[n + 1] = Cc_star_col_map[0]-1;
      else
        Cc_star_col_map[n + 1] = Cc_star_col_map[0]+1;

    std::vector<int64_t> Cc_star_col_map4knowns(Nn + 1, 0);

    Cc_star_col_map4knowns[0] = static_cast<int64_t>(c);
    for (size_t n=0; n<Nn; ++n)
    {
      const int64_t nmap = pwl.MapDOFLocal(cell_c, n) + static_cast<int64_t>(Nc);
      Cc_star_col_map4knowns[n + 1] = nmap;
    }

    //TODO: End : Find generalization of mapping

    //=========================================== Assemble connectivity
    const double D_c_n   = -speed_of_light_cmpsh / (3 * sigma_t_c_n);
    const double D_c_np1 = -speed_of_light_cmpsh / (3 * sigma_t_c_np1);
    double grad_dot_J_n = 0.0;
    //====================================== Loop over faces
    const auto& Cc_star = list_of_Cc_star[c];
    for (size_t f=0; f<Nf; ++f)
    {
      const auto&  face = cell_c.faces[f];
      const auto&  nf   = face.normal;
      const auto&  Sf   = S_surf[f];
      const size_t Nfn  = face.vertex_ids.size();

      //=============================== Assemble primary equation portion
      const auto& col_map       = Cc_star_col_map;
      const auto& known_col_map = Cc_star_col_map4knowns;
      for (size_t fj=0; fj<Nfn; ++fj)
      {
        const size_t  j    = fem_data.FaceDofMapping(f, fj);

        for (size_t d = 0; d < Nd; ++d)
          for (size_t n = 0; n < (Nn + 1); ++n)
          {
            A[col_map[0]][col_map[n]] +=
              nf[d] * Sf[j] * (theta1*D_c_np1/V_c) * Cc_star[j * Nd + d][n];
            grad_dot_J_n +=
              nf[d] * Sf[j] * (D_c_n/V_c) * Cc_star[j * Nd + d][n] *
              rad_E_n[known_col_map[n]];
          }
      }

      //=============================== Assembly auxiliary portions
      for (size_t fj=0; fj<Nfn; ++fj)
      {
        const size_t  j    = fem_data.FaceDofMapping(f, fj);

        for (size_t d=0; d<Nd; ++d)
          for (size_t n=0; n<(Nn+1); ++n)
            A[col_map[j+1]][col_map[n]] -=
              nf[d] * Sf[j] * (theta1*D_c_np1/V_c) * Cc_star[j * Nd + d][n];

      }//for fj
    }//for f

    //=========================================== Diagonal and rhs
    const auto& col_map = Cc_star_col_map;

    A[col_map[0]][col_map[0]] += tau + k1 + k4*k5;
    b[col_map[0]]             += -k3 - k4*k6 + tau*rad_E_c_nphstar
                                 - theta2*grad_dot_J_n;
  }//for cell

  //============================================= Solve system and reconstruct
  auto x = chi_math::BandedSolver(A,b,2,2);

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id;
    rad_E_np1[c] = x[2*c+1];

    for (int j=0; j<cell.vertex_ids.size(); ++j)
    {
      const int64_t jmap = pwl.MapDOFLocal(cell, j) + static_cast<int64_t>(Nc);

      if (j==0) rad_E_np1[jmap] = x[2*c];
      else      rad_E_np1[jmap] = x[2*c+2];
    }
  }//for cell

  InverseAlgebraFor_E(sim_refs, rad_E_np1, U_np1);
}