#include "radtranGreyDiff.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

void chi_radhydro::RadTranGreyDiffusion::Execute()
{
  chi::log.Log0() << "\nExecuting " << this->TextName() << " solver\n\n";

  /**Lambda to map e in the e-rad_E system.*/
  auto Map_e = [this](const chi_mesh::Cell& cell)
  {
    return fv->MapDOFLocal(cell,0,e_E_uk_man,0,0);
  };

  /**Lambda to map rad_E in the e-rad_E system.*/
  auto Map_rad_E = [this](const chi_mesh::Cell& cell)
  {
    return fv->MapDOFLocal(cell,0,e_E_uk_man,1,0);
  };

  std::vector<UVector>&             U_n = U_old;
  std::vector<UVector>              U_n_star(num_nodes_local);
  std::vector<UVector>              U_nph(num_nodes_local);
  std::vector<UVector>              U_nph_star(num_nodes_local);
  std::vector<UVector>              U_np1(num_nodes_local);

  std::vector<double>&              rad_E_n = rad_E;
  std::vector<double>               rad_E_n_star(num_nodes_local);
  std::vector<double>               rad_E_nph(num_nodes_local);
  std::vector<double>               rad_E_nph_star(num_nodes_local);
  std::vector<double>               rad_E_np1(num_nodes_local);

  FieldsToU(U_n, *grid, rho,u,v,w,e);

  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double delta_t_hydro =
      ComputeCourantLimitDelta_t(*grid, rho,u,v,w,p,gamma,cell_char_length,C_cfl);
    const double dt            = std::min(delta_t_max, delta_t_hydro);

    chi::log.Log() << "Timestep " << n << " with dt=" << dt << " t=" << time;

    //##################################################### PREDICTOR
    //================================= Compute kappas
    VecDbl kappa_t(num_nodes_local,0.0);
    VecDbl kappa_a(num_nodes_local,0.0);

    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c = cell.local_id;

      const double kappa_s = 0.0;

      const double kappa_1 = 577.35;
      const double kappa_2 = 0.0;
      const double kappa_3 = 1.0;
      const double n_exponent = 0.0;
      const double T = temperature[c];
      kappa_a[c] = kappa_1/(kappa_2*pow(T,n_exponent) + kappa_3);

      kappa_t[c] = kappa_s + kappa_a[c];
    }//for cell

    //================================= Compute gradients
    //Populates grad_U with slope limiters
    const auto grad_U     = ComputeUGradients(U_n);
    const auto grad_rad_E = ComputeRadEGradients(rad_E_n);

    //================================= Advance half timestep
    // Advection of U and rad_E
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;
      const Vec3&    x_cc    = cell.centroid;

      const UVector U_c_n    = U_n[c];
      const double  rad_E_c_n = rad_E_n[c];

      UVector U_c_n_star     = U_c_n;
      double  rad_E_c_n_star = rad_E_c_n;

      const double   tau     = 0.5*dt;

      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f  = fv_view->face_area[f];
        const Vec3&  n_f  = cell.faces[f].normal;
        const Vec3&  x_fc = cell.faces[f].centroid;

        const UVector  U_f   = UplusDXGradU(U_c_n, x_fc-x_cc, grad_U[c]);
        const Vec3&    u_f   = Vec3(U_f[1],U_f[2],U_f[3])/U_f[0];
        const double rad_E_f = rad_E_c_n + (x_fc-x_cc).Dot(grad_rad_E[c]);

        const UVector F_f = MakeF(U_f, p[c], n_f);

        U_c_n_star     -= (1/tau) * (1/V_c) * A_f * F_f;
        rad_E_c_n_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_E_f * u_f);
      }//for f

      U_n_star[c]     = U_c_n_star;
      rad_E_n_star[c] = rad_E_c_n_star;
    }//for cell

    // Update density and momentum to nph
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;

      const double   tau     = 0.5*dt;

      const double D_c     = -speed_of_light_cmpsh / (3 * rho[c] * kappa_t[c]);
      const double rad_E_c = rad_E_n[c];

      UVector U_c_nph_01 = U_n_star[c];

      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const Vec3 A_f  = fv_view->face_area[f] * cell.faces[f].normal;
        const Vec3 x_cf = cell.faces[f].centroid - cell.centroid;

        const double k_c = D_c/x_cf.Norm();

        double k_cn = k_c;
        double rad_E_cn = rad_E_c;
        if (cell.faces[f].has_neighbor)
        {
          const uint64_t cn = cell.faces[f].neighbor_id;
          const auto& adj_cell = grid->cells[cn];
          const Vec3   x_fcn = adj_cell.centroid - cell.faces[f].centroid;

          const double D_cn = -speed_of_light_cmpsh / (3 * rho[cn] * kappa_t[cn]);

          k_cn = D_cn/x_fcn.Norm();
          rad_E_cn = rad_E_n[cn];
        }

        const double rad_E_f = (k_cn*rad_E_cn + k_c*rad_E_c)/
                               (k_cn + k_c);
        U_c_nph_01(1) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.x * rad_E_f;
        U_c_nph_01(2) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.y * rad_E_f;
        U_c_nph_01(3) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.z * rad_E_f;
      }//for f

      U_nph[c] = U_c_nph_01;
    }//for cell

    // Internal energy and radiation energy
    {
      const auto num_e_E_dofs = fv->GetNumLocalDOFs(e_E_uk_man);
      MatDbl A(num_e_E_dofs, VecDbl(num_e_E_dofs,0.0));
      VecDbl b(num_e_E_dofs, 0.0);

      for (const auto& cell_c : grid->local_cells)
      {
        const uint64_t c = cell_c.local_id;
        const double T_c_n      = IdealGasTemperatureFromCellU(U_n[c], Cv[c]);
        const double T_c_n_star = IdealGasTemperatureFromCellU(U_n_star[c], Cv[c]);
        const double rho_c_nph  = U_nph[c][0];
        const double sigma_tc_nph = rho_c_nph * kappa_t[c];
        const double sigma_ac_nph = rho_c_nph * kappa_a[c];

        const uint64_t cmap_e = Map_e(cell_c);
        const uint64_t cmap_E = Map_rad_E(cell_c);


      }//for cell
    }

    //##################################################### CORRECTOR
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;
      const Vec3&    x_cc    = cell.centroid;

      const UVector& U_c_n   = U_n[c];
      UVector  U_c_np1 = U_c_n;

      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f  = fv_view->face_area[f];
        const Vec3&  n_f  = cell.faces[f].normal;
        const Vec3&  x_fc = cell.faces[f].centroid;

        const UVector& U_L     = UplusDXGradU(U_nph[c], x_fc-x_cc, grad_U[c]);
        const double   gamma_L = gamma[c];
        const double   p_L     = IdealGasPressureFromCellU(U_L, gamma_L);

        UVector  U_R     = U_L;
        double   gamma_R = gamma_L;
        double   p_R     = p_L;

        if (cell.faces[f].has_neighbor)
        {
          const uint64_t cn            = cell.faces[f].neighbor_id;
          const auto&    adjacent_cell = grid->cells[cn];
          const Vec3&    adj_x_cc       = adjacent_cell.centroid;

          U_R     = UplusDXGradU(U_nph[cn], x_fc-adj_x_cc, grad_U[cn]);
          gamma_R = gamma[cn];
          p_R     = IdealGasPressureFromCellU(U_R, gamma_R);
        }

        const FVector F_hllc_f = HLLC_RiemannSolve(U_L    , U_R,
                                                   p_L    , p_R,
                                                   gamma_L, gamma_R,
                                                   n_f);

        U_c_np1 -= (dt/V_c)* A_f * F_hllc_f;
      }//for f
      U_np1[c] = U_c_np1;
    }//for cell

    U_n = U_np1;
    UToFields(U_n,*grid,rho,u,v,w,e,p,gamma);

    time += dt;
    if (t_max >= 0.0 and time >= t_max)
      break;
  }//for n timesteps

  chi::log.Log0() << "\nDone executing " << this->TextName() << " solver\n\n";
}

