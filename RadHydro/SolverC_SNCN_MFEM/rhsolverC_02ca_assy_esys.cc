#include "rhsolverC.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "ChiMath/chi_math_banded_solvers.h"

void chi_radhydro::SolverC_SNCN_MFEM::
  AssySolveMFEMEnergyExchange(
  SimRefs&                        sim_refs,
  SDM_PWLC&                       pwl,
  SDM_PWLD&                       pwld,
  std::vector<MatDbl>&            list_of_Cc_star,
  std::vector<MatDbl>&            list_of_Cvol_star,
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
  const std::vector<Vec3>&        grad_rad_E_nph,
  const std::vector<Vec3>&        rad_F0_nph,
  const std::vector<Vec3>&        rad_F_n,
  const std::vector<double>&      VEFf_nodal,
  const std::vector<double>&      VEFf_ctr,
  std::vector<UVector>&           U_np1,
  std::vector<double>&            rad_E_np1,
  std::vector<Vec3>&              rad_F_np1,
  std::vector<Vec3>&              rad_F0_np1)
{
//  //============================================= Define lambda's
//  auto PrintSystem = [](const MatDbl& mat, const VecDbl& rhs={},
//                        const double div_fac_mat=1.0,
//                        const double div_fac_rhs=1.0)
//  {
//    int j=0;
//    for (const auto& row : mat)
//    {
//      for (double val : row)
//        if (std::fabs(val) > 1.0e-12)
//          printf("%8.1e ",val/div_fac_mat);
//        else
//          printf("00000000 ");
//      if (not rhs.empty())
//        printf("     %5.1f ",rhs[j]/div_fac_rhs);
//      std::cout << "\n";
//      ++j;
//    }
//  };

  const auto&  grid = sim_refs.grid;
        auto&  fv   = sim_refs.fv;
  const double Cv   = sim_refs.Cv;

  const size_t Nc    = grid.local_cells.size(); //Number of cells
  const size_t Ns    = rad_E_n.size();
  const double c_spd = chi_radhydro::speed_of_light_cmpsh;

  const size_t Nd = 3; //Number of dimensions

  MatDbl A;
  VecDbl b;
  A.assign(Ns, VecDbl(Ns, 0.0));
  b.assign(Ns, 0.0);

  sim_refs.e_recon_coeff_k5.assign(Ns, 0.0);
  sim_refs.e_recon_coeff_k6.assign(Ns, 0.0);

  std::vector<std::array<double,4>> list_of_a_coeffs(Nc, {0.0,0.0,0.0,0.0});

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
    const double   rho_c_nph   = U_nph[c][RHO];
    const double   rho_c_np1   = U_np1[c][RHO];

    const Vec3     u_c_nph     = VelocityFromCellU(U_nph[c]);

    const double   T_c_n       = IdealGasTemperatureFromCellU(U_n[c], Cv);
    const double   T_c_nphstar = IdealGasTemperatureFromCellU(U_nphstar[c], Cv);

    const double   E_c_nphstar = U_nphstar[c][MAT_E];

    const double   e_c_nphstar     = InternalEnergyFromCellU(U_nphstar[c]);
    const double   rad_E_c_n       = rad_E_n[c];
    const double   rad_E_c_nph     = rad_E_nph[c];

    const Vec3     u_c_np1         = VelocityFromCellU(U_np1[c]);
    const double   VEFf_c           = VEFf_ctr[c];

    //=========================================== Compute sigmas
    const double sigma_t_c_n   = rho_c_n * kappa_t_n[c];
    const double sigma_a_c_n   = rho_c_n * kappa_a_n[c];

    const double sigma_t_c_nph = rho_c_nph * kappa_t_nph[c];

    const double sigma_t_c_np1 = rho_c_np1 * kappa_t_nph[c];
    const double sigma_a_c_np1 = rho_c_np1 * kappa_a_nph[c];

    //=========================================== Catch infinite diffusion coeffs
    if (std::fabs(sigma_t_c_np1) < 1.0e-6)
    {
      A[c][c] = 1.0; b[c] = rad_E_c_nph;
      sim_refs.e_recon_coeff_k5[c] = 0.0;
      sim_refs.e_recon_coeff_k6[c] = e_c_nphstar;
      continue;
    }

    //=========================================== Compute explicit terms
    Vec3 rad_F0_c_nph;
    for (size_t i=0; i<Nn; ++i)
    {
      const auto imap = pwld.MapDOFLocal(cell_c, i);
      rad_F0_c_nph += rad_F0_nph[imap];
    }
    rad_F0_c_nph /= static_cast<double>(Nn);
    double Srp_nph = (sigma_t_c_nph / c_spd) * rad_F0_c_nph.Dot(u_c_nph);

    // sigma_a c (aT^4 - radE)
    double Sea_n = MakeEmAbsSource(sigma_a_c_n, T_c_n, rad_E_c_n);

    //=========================================== Compute constants
    const double k1 = theta1 * sigma_a_c_np1 * c_spd;
    const double k2 = 4.0 * pow(T_c_nphstar, 3.0) / Cv;

    const double k3 = - k1*a*pow(T_c_nphstar, 4.0)
                      + k1*a*k2*e_c_nphstar
                      - theta2 * Sea_n
                      + Srp_nph;
    const double k4 = -k1*a*k2;

    const double k5 = k1/(tau * rho_c_np1 - k4);
    const double k6 =
      (k3 - tau * 0.5 * rho_c_np1 * u_c_np1.NormSquare() + tau * E_c_nphstar) /
      (tau * rho_c_np1 - k4);

    sim_refs.e_recon_coeff_k5[c] = k5;
    sim_refs.e_recon_coeff_k6[c] = k6;

    //====================================== Nodal map the columns of Cc_star
    const auto& Cc_star   = list_of_Cc_star[c];
    const auto& Cvol_star = list_of_Cvol_star[c];
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
    const double c_sqrd = pow(c_spd, 2.0);
    const double a1 = theta1 * c_sqrd / tau;
    const double a2 = theta2 * c_sqrd / tau;
    const double a3 = sigma_t_c_nph * c_spd / tau;

    const double divisor = (1.0 + (a1/c_spd)*sigma_t_c_np1);
    const double a4 = a1/divisor;
    const double a5 = (1.0 - (a2/c_spd)*sigma_t_c_n)/divisor;
    const double a6 = a2/divisor;
    const double a7 = a3/divisor;

    list_of_a_coeffs[c] = {a4,a5,a6,a7};


    double grad_dot_F_n = 0.0;
    //====================================== Loop over faces
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
        const int64_t jmap = pwld.MapDOFLocal(cell_c, j);

        //Implicit fE term
        for (size_t d=0; d<Nd; ++d)
          for (size_t n=0; n <(Nn + 1); ++n)
          {
            const double VEFf_n = (n == 0)?
                                  VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                  VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
            A[col_map[0]][col_map[n]] +=
              (theta1 * (-a4)/V_c) *
              nf[d] * Sf[j] * Cc_star[j * Nd + d][n] *
              VEFf_n;
          }//for n

        //Explicit F term
        for (size_t d=0; d<Nd; ++d)
          b[col_map[0]] -= (theta1 * a5/V_c)*nf[d]*Sf[j]*rad_F_n[jmap][d];

        //Explicit fE term
        for (size_t d=0; d<Nd; ++d)
          for (size_t n=0; n <(Nn + 1); ++n)
          {
            const double VEFf_n = (n == 0)?
                                  VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                  VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
            b[col_map[0]] -=
              (theta1 * (-a6)/V_c) *
              nf[d] * Sf[j] * Cc_star[j * Nd + d][n] *
              VEFf_n * rad_E_n[known_col_map[n]];
          }

        //Explicit convection term
        for (size_t d=0; d<Nd; ++d)
          b[col_map[0]] -=
            (theta1 * a7/V_c) * nf[d] * Sf[j] * Cvol_star[j * Nd + d][d] *
            (1.0 + VEFf_c) * rad_E_c_nph * u_c_nph[d];

        for (size_t d=0; d<Nd; ++d)
          grad_dot_F_n += (1.0/V_c) * nf[d] * Sf[j] * rad_F_n[jmap][d];
      }

      //=============================== Assemble auxiliary portions
      for (size_t fj=0; fj<Nfn; ++fj)
      {
        const size_t  j    = fem_data.FaceDofMapping(f, fj);
        const int64_t jmap = pwld.MapDOFLocal(cell_c, j);

        //Implicit fE term
        for (size_t d=0; d<Nd; ++d)
          for (size_t n=0; n <(Nn + 1); ++n)
          {
            const double VEFf_n = (n == 0)?
                                  VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                  VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
            A[col_map[j+1]][col_map[n]] +=
              (-a4) *
              nf[d] * Sf[j] * Cc_star[j * Nd + d][n] *
              VEFf_n;
          }

        //Explicit F term
        for (size_t d=0; d<Nd; ++d)
          b[col_map[j+1]] += -a5*nf[d]*Sf[j]*rad_F_n[jmap][d];

        //Explicit fE term
        for (size_t d=0; d<Nd; ++d)
          for (size_t n=0; n <(Nn + 1); ++n)
          {
            const double VEFf_n = (n == 0)?
                                  VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                  VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
            b[col_map[j+1]] +=
              a6 *
              nf[d] * Sf[j] * Cc_star[j * Nd + d][n] *
              VEFf_n * rad_E_n[known_col_map[n]];
          }

        //Explicit convection term
        for (size_t d=0; d<Nd; ++d)
          b[col_map[j+1]] +=
            (-a7) * nf[d] * Sf[j] * Cvol_star[j * Nd + d][d] *
            (1.0 + VEFf_c) * rad_E_c_nph * u_c_nph[d];

      }//for fj

      if (not face.has_neighbor)
      {
        const double radE_adj = chi_radhydro::MakeRadEFromBC(
          sim_refs.bc_settings.at(face.neighbor_id), rad_E_n[c]);
        const UVector U_adj = chi_radhydro::MakeUFromBC(
          sim_refs.bc_settings.at(face.neighbor_id), U_n[c]);
        const Vec3 u_adj = VelocityFromCellU(U_adj);

        for (size_t d=0; d<Nd; ++d)
          b[col_map[f+1]] += (4.0/3.0) * radE_adj * u_adj[d]/c_spd;
      }
    }//for f

    //=========================================== Diagonal and rhs
    const auto& col_map = Cc_star_col_map;

    A[col_map[0]][col_map[0]] += tau + k1 + k4*k5;
    b[col_map[0]]             += -k3 - k4*k6 + tau*rad_E_c_n
                                 - theta2 * grad_dot_F_n;
  }//for cell

  //============================================= Solve system
  auto x = chi_math::BandedSolver(A,b,2,2);

  //============================================= Reconstruct RadE
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

  //============================================= Compute RadF and RadF0
  for (const auto& cell_c : grid.local_cells)
  {
    //================================= Get cell geometry info
    const uint64_t c         = cell_c.local_id;

    //================================= Get fem stuff
    const auto&    fem_data = pwl.GetUnitIntegrals(cell_c);
    const size_t   Nn       = fem_data.NumNodes();

    const auto&    S_surf   = fem_data.GetIntS_shapeI();

    //================================= Get cell physics-info
    const double VEFf_c  = VEFf_ctr[c];
    const Vec3   u_c_nph = VelocityFromCellU(U_nph[c]);

    //================================= Get construction coefficients
    const auto& a_coeffs = list_of_a_coeffs[c];
    const double a4 = a_coeffs[0];
    const double a5 = a_coeffs[1];
    const double a6 = a_coeffs[2];
    const double a7 = a_coeffs[3];

    //================================= Map knowns
    const auto& Cc_star = list_of_Cc_star[c];
    const auto& Cvol_star = list_of_Cvol_star[c];
    std::vector<int64_t> Cc_star_col_map4knowns(Nn + 1, 0);

    Cc_star_col_map4knowns[0] = static_cast<int64_t>(c);
    for (size_t n=0; n<Nn; ++n)
    {
      const int64_t nmap = pwl.MapDOFLocal(cell_c, n) + static_cast<int64_t>(Nc);
      Cc_star_col_map4knowns[n + 1] = nmap;
    }

    //================================= Compute F0 and F
    const auto& known_col_map = Cc_star_col_map4knowns;
    for (size_t i=0; i<Nn; ++i)
    {
      const int64_t imap = pwld.MapDOFLocal(cell_c, i);
      auto& radF_c_i_np1  = rad_F_np1[imap];
      auto& radF0_c_i_np1 = rad_F0_np1[imap];

      radF_c_i_np1 = {0.0,0.0,0.0};
      radF0_c_i_np1 = {0.0,0.0,0.0};

      // fE^n+1
      for (size_t d=0; d<Nd; ++d)
        for (size_t n=0; n<(Nn+1); ++n)
        {
          const double VEFf_n = (n == 0)?
                                VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
          radF0_c_i_np1(d) += -a4 * Cc_star[i*Nd + d][n] *
                              VEFf_n * rad_E_np1[known_col_map[n]];
        }

      // F^n
      radF0_c_i_np1 += a5 * rad_F_n[imap];

      // fE^n
      for (size_t d=0; d<Nd; ++d)
        for (size_t n=0; n<(Nn+1); ++n)
        {
          const double VEFf_n = (n == 0)?
                                VEFf_ctr  [Cc_star_col_map4knowns[n]] :
                                VEFf_nodal[Cc_star_col_map4knowns[n]-Nc];
          radF0_c_i_np1(d) += -a6 * Cc_star[i*Nd + d][n] *
                              VEFf_n * rad_E_n[known_col_map[n]];
        }

      radF_c_i_np1 = radF0_c_i_np1;

      // Convection term
      for (size_t d=0; d<Nd; ++d)
        radF_c_i_np1(d) +=
          a7 * Cvol_star[i*Nd + d][d] *
          (1.0 + VEFf_c) * rad_E_nph[c] * u_c_nph[d];

    }//for i
  }//for cell

  InverseAlgebraFor_E(sim_refs, rad_E_np1, U_np1);
}