#include "rhsolverC.h"

#include "ChiMath/chi_math_range.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

void chi_radhydro::SolverC_SNCN_MFEM::
  Sweep1D(SimRefs&                           sim_refs,
          SDM_PWLC&                          pwlc,
          SDM_PWLD&                          pwld,
          const std::vector<double>&         kappa_a,
          const std::vector<double>&         kappa_t,
          const std::vector<UVector>&        U,
          const std::vector<double>&         radE,
          const std::vector<Vec3>&           radF0,
          const chi_math::AngularQuadrature& quadrature,
          VecDbl& VEFf_nodal,
          VecDbl& VEFf_ctr,
          const bool verbose/*=false*/)
{
  const std::string fname = __FUNCTION__;
  const double ONE_DIV_4PI = 1.0/(4.0*M_PI);
  const double c_spd = speed_of_light_cmpsh;

  //============================================= Lambda
  auto GetUpstreamPsi = [&sim_refs, &pwld, &c_spd, &ONE_DIV_4PI](
    const std::vector<double>& psi,
    const chi_mesh::CellFace& face,
    const size_t face_index,
    const Vec3& omega,
    const double T_to_4th,
    const Vec3& u)
  {
    const uint64_t neighbor_id = face.neighbor_id;

    if (face.has_neighbor)
    {
      const auto& neighbor_cell = sim_refs.grid.local_cells[neighbor_id];

      int neighbor_node;
      if (face_index == 0) neighbor_node = 1;
      if (face_index == 1) neighbor_node = 0;
      const int64_t node_map = pwld.MapDOFLocal(neighbor_cell, neighbor_node);

      return psi[node_map];
    }
    else
    {
      return (1.0 + 4.0*omega.Dot(u)/c_spd) *
             a * c_spd * T_to_4th * ONE_DIV_4PI;
    }
  };

  //============================================= Grid check
  const auto& grid = sim_refs.grid;

  if (not (grid.Attributes() & chi_mesh::DIMENSION_1))
    throw std::logic_error(fname + ": Can only be used for 1D meshes.");

  if (not (grid.Attributes() & chi_mesh::ORTHOGONAL))
    throw std::logic_error(fname + ": Can only be used for 1D "
                                   "orthogonal meshes.");

  //============================================= Initializing of Misc items
  const auto ONE_DOF_PER_NODE =
    chi_math::UnknownManager::GetUnitaryUnknownManager();

  const size_t Nc = grid.local_cells.size();
  const int Nc_int = static_cast<int>(Nc);
  const size_t num_pwld_nodes = pwld.GetNumLocalDOFs(ONE_DOF_PER_NODE);
  const size_t num_pwlc_nodes = pwlc.GetNumLocalDOFs(ONE_DOF_PER_NODE);

  MatDbl Amat (2, VecDbl(2,0.0));
  MatDbl Atemp(2, VecDbl(2,0.0));
  VecDbl b    (2, 0.0);

  VecDbl psi(num_pwld_nodes, 0.0);
  VecDbl phi(num_pwld_nodes, 0.0);
  VecDbl VEFftop_nodal(num_pwlc_nodes, 0.0);
  VecDbl VEFfbot_nodal(num_pwlc_nodes, 0.0);
  VecDbl VEFftop_ctr(Nc, 0.0);
  VecDbl VEFfbot_ctr(Nc, 0.0);

  if (verbose)
    chi::log.Log() << "Done preparing for sweeps.";

  //============================================= Loop over angles
  const size_t ND = quadrature.omegas.size();
  for (size_t n=0; n<ND; ++n)
  {
    const auto& omega  = quadrature.omegas[n];
    const auto& weight = quadrature.weights[n];
    const double mu_omega = omega.z;
    const double mu_sqr = mu_omega * mu_omega;

    if (verbose)
      chi::log.Log() << "Sweeping direction " << n << " " << omega.PrintStr();

    //====================================== Determine sweep order
    std::vector<int> sweep_order;
    if (omega.z > 0.0) sweep_order = chi_math::Range<int>(0, Nc_int);
    else               sweep_order = chi_math::Range<int>(Nc_int-1,-1,-1);


    //====================================== Loop over cells
    psi.assign(num_pwld_nodes, 0.0);
    for (const int cell_local_id : sweep_order)
    {
      //=============================== Get geometry info
      const auto& cell = grid.local_cells[cell_local_id];
      const uint64_t c = cell.local_id;

      const auto& fe_intgrl_vals = pwld.GetUnitIntegrals(cell);

      const size_t Nn = fe_intgrl_vals.NumNodes();
      const size_t Nf = cell.faces.size();

      const auto& G      = fe_intgrl_vals.GetIntV_shapeI_gradshapeJ();
      const auto& M      = fe_intgrl_vals.GetIntV_shapeI_shapeJ();
      const auto& M_surf = fe_intgrl_vals.GetIntS_shapeI_shapeJ();
      const auto& V_ci   = fe_intgrl_vals.GetIntV_shapeI();

      //=============================== Get physics quantities
      const double T = IdealGasTemperatureFromCellU(U[c], sim_refs.Cv);
      const double T_to_4th = pow(T,4.0);
      const double rho = U[c][RHO];

      const double sigma_s = rho * (kappa_t[c] - kappa_a[c]);
      const double sigma_a = rho * kappa_a[c];
      const double sigma_t = sigma_a + sigma_s;
      const double radE_c = radE[c];
      const Vec3   u = VelocityFromCellU(U[c]);

      Amat.assign(2, VecDbl(2,0.0));
      b.assign(2,0.0);
      if (verbose)
        chi::log.Log() << "sigma_t " << sigma_t << " " << 1/(V_ci[0] + V_ci[1]);

      //=============================== Surface integrals
      for (size_t f=0; f<Nf; ++f)
      {
        const auto& face = cell.faces[f];
        const double mu = face.normal.Dot(omega);
        if (mu < 0.0)
        {
          const double us_psi = GetUpstreamPsi(psi, face, f, omega, T_to_4th, u);
          const size_t i=f, j=f;

          const double mu_Mf_ij = -mu * M_surf[f][i][j];
          Amat[i][j] += mu_Mf_ij;
          b[i]       += mu_Mf_ij * us_psi;

          if (verbose)
            chi::log.Log() << "us_psi[" << c << "]=" << us_psi;
        }//if incident
      }//for f

      if (verbose)
        chi::log.Log() << "surface_B[" << c << "] " << b[0] << " " << b[1];

      //=============================== Complete source term
      // q = M_n^T * q_moms
      VecDbl src_b(2,0.0);
      for (size_t i=0; i<Nn; ++i)
      {
        double b_i = V_ci[i] * (  ONE_DIV_4PI * sigma_s * c_spd * radE_c
                                + ONE_DIV_4PI * sigma_a * a * c_spd * T_to_4th
                                + (sigma_t/M_PI) * radE_c * omega.Dot(u));
        for (size_t j=0; j<Nn; ++j)
        {
          const int64_t jmap = pwld.MapDOFLocal(cell, j);
          const auto& radF0_j_dot_u = radF0[jmap].Dot(u);

          b_i -= ONE_DIV_4PI * (sigma_t / c_spd) * M[i][j] * radF0_j_dot_u;
        }

        b[i] += b_i;
        src_b[i] = b_i;
      }//for i

      if (verbose)
        chi::log.Log() << "src_B    [" << c << "] " << src_b[0] << " " << src_b[1];

      //=============================== Add G and M to Amat
      for (int i = 0; i < Nn; ++i)
        for (int j = 0; j < Nn; ++j)
          Amat[i][j] += omega.Dot(G[i][j]) + M[i][j] * sigma_t;

      //=============================== Solve
      chi_math::GaussElimination(Amat, b, 2);
      const auto& cell_psi = b;

      if (verbose)
        chi::log.Log() << "cell_psi [" << c << "] " << cell_psi[0] << " " << cell_psi[1];

      //=============================== Store psi and contribute to phi
      for (size_t i=0; i<Nn; ++i)
      {
        const int64_t imap = pwld.MapDOFLocal(cell, i);
        psi[imap]  = cell_psi[i];
        phi[imap] += weight * cell_psi[i];
      }//for i

      //=============================== Compute VEF-factors
      for (size_t f=0; f<Nf; ++f)
      {
        const auto& face = cell.faces[f];
        const int64_t imap = pwlc.MapDOFLocal(cell,f);

        if (omega.Dot(face.normal) > 0.0)
        {
          VEFftop_nodal[imap] += weight * mu_sqr * cell_psi[f];
          VEFfbot_nodal[imap] += weight * cell_psi[f];
        }
      }//for f

      VEFftop_ctr[c] += weight * mu_sqr * 0.5 * (cell_psi[0] + cell_psi[1]);
      VEFfbot_ctr[c] += weight *          0.5 * (cell_psi[0] + cell_psi[1]);
    }//for cell in sweep order
  }//for omega_n

  if (verbose)
    chi::log.Log() << "Done sweeping.";

  //============================================= Finish computing VEFf
  for (size_t n=0; n<num_pwlc_nodes; ++n)
    VEFftop_nodal[n] /= VEFfbot_nodal[n];

  for (size_t c=0; c<Nc; ++c)
    VEFftop_ctr[c] /= VEFfbot_ctr[c];

  VEFf_nodal = std::move(VEFftop_nodal);
  VEFf_ctr   = std::move(VEFftop_ctr);

}