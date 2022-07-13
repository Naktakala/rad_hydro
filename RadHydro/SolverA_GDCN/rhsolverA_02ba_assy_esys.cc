#include "rhsolverA.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

void chi_radhydro::SolverA_GDCN::
  AssembleGeneralEnergySystem(
  const chi_mesh::MeshContinuum&  grid_ref,
  std::shared_ptr<SDM_FV>&        fv_ref,
  const std::map<uint64_t, BCSetting>& bc_setttings,
  const std::vector<double>&      kappa_a_n,
  const std::vector<double>&      kappa_t_n,
  const std::vector<double>&      kappa_a_nph,
  const std::vector<double>&      kappa_t_nph,
  const double                    Cv,
  double                          tau,
  double                          theta1,
  double                          theta2,
  const std::vector<UVector>&     U_n,
  const std::vector<UVector>&     U_nph,
  const std::vector<UVector>&     U_nphstar,
  const std::vector<UVector>&     U_np1,
  const std::vector<GradUTensor>& grad_U_nph,
  const std::vector<double>&      rad_E_n,
  const std::vector<double>&      rad_E_nph,
  const std::vector<double>&      rad_E_nphstar,
  const std::vector<Vec3>&        grad_rad_E_nph,
  std::vector<double>&            k5_vec,
  std::vector<double>&            k6_vec,
  MatDbl &A, VecDbl &b)
{
  const auto& grid    = grid_ref;
  const auto& fv      = fv_ref;

  k5_vec.assign(grid.local_cells.size(), 0.0);
  k6_vec.assign(grid.local_cells.size(), 0.0);

  for (const auto& cell_c : grid.local_cells)
  {
    //=========================================== Get cell geometry info
    const uint64_t c         = cell_c.local_id;
    const size_t   num_faces = cell_c.faces.size();

    const auto&    fv_view_c   = fv->MapFeView(c);
    const auto&    face_areas  = fv_view_c->face_area;
    const Vec3&    x_c         = cell_c.centroid;
    const double   V_c         = fv_view_c->volume;

    //=========================================== Get cell physics-info
    const double   rho_c_n     = U_n[c][RHO];
    const double   rho_c_np1   = U_np1[c][RHO];
    const double   T_c_n       = IdealGasTemperatureFromCellU(U_n[c], Cv);
    const double   T_c_nph     = IdealGasTemperatureFromCellU(U_nph[c], Cv);
    const double   T_c_nphstar = IdealGasTemperatureFromCellU(U_nphstar[c], Cv);

    const double   E_c_nphstar = U_nphstar[c][MAT_E];

    const double   e_c_nphstar     = InternalEnergyFromCellU(U_nphstar[c]);
    const double   rad_E_c_n       = rad_E_n[c];
    const double   rad_E_c_nph     = rad_E_nph[c];
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
      k5_vec[c] = 0.0; k6_vec[c] = e_c_nphstar;
      continue;
    }

    //=========================================== Compute explicit terms
    double third_grad_rad_E_dot_u_nph =
      Make3rdGradRadE(cell_c, V_c, face_areas,
                      rad_E_c_nph, grad_rad_E_nph[c],
                      U_nph[c], grad_U_nph[c]);

    // sigma_a c (aT^4 - radE)
    double Sea_n = MakeEmAbsSource(sigma_a_c_n, T_c_n, rad_E_c_n);

    double grad_dot_J_n = ComputeGradDotJ(grid, fv_ref, cell_c,
                                          sigma_t_c_n, kappa_t_n,
                                          U_n, rad_E_n);

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

    k5_vec[c] = k5;
    k6_vec[c] = k6;

    //=========================================== Assemble connectivity
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell_c.faces[f];
      const auto& x_f  = face.centroid;
      const auto  A_f  = face_areas[f] * face.normal;

      double sigma_t_cn_np1 = sigma_t_c_np1;
      Vec3   x_cn          = x_c + 2*(x_f-x_c);

      if (not face.has_neighbor) //DEFAULT REFLECTING BC for radE
      {
        //Nothing to do for reflecting bc
      }
      else                       //NEIGHBOR CELL
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid.cells[cn];
        x_cn = cell_cn.centroid;

        const double rho_cn_np1 = U_np1[cn][0];

        sigma_t_cn_np1 = rho_cn_np1 * kappa_t_nph[cn];
      }

      const double sigma_tf_np1 = (sigma_t_c_np1 + sigma_t_cn_np1) / 2.0;

      const double Df_np1 = -speed_of_light_cmpsh / (3.0 * sigma_tf_np1);

      const Vec3 x_ccn = x_cn - x_c;
      const Vec3 kf_np1 = Df_np1 * x_ccn / x_ccn.NormSquare();

      const double coeff_LHS = (theta1 / V_c) * A_f.Dot(kf_np1);

      if (not face.has_neighbor) //DEFAULT REFLECTING
      {
        //J_f = 0 therefor no connectivity elements
      }
      else                       //NEIGHBOR CELL
      {
        A[c][c] += -coeff_LHS;
        const uint64_t cn = face.neighbor_id;
        A[c][cn] += coeff_LHS;
      }
    }//for f in connectivity

//    for (size_t f=0; f<num_faces; ++f)
//    {
//      const auto& face = cell_c.faces[f];
//      const auto& x_f  = face.centroid;
//      const auto  A_f  = face_areas[f] * face.normal;
//      const auto  x_cf = x_f - x_c;
//
//      const double D_c = -speed_of_light_cmpsh / (3 * sigma_t_c_np1);
//      const double k_c = D_c/x_cf.Norm();
//
//      double sigma_t_cn_np1 = sigma_t_c_np1;
//      Vec3   x_cn          = x_c + 2*(x_f-x_c);
//
//      if (not face.has_neighbor) //DEFAULT REFLECTING BC for radE
//      {
//        //Nothing to do for reflecting bc
//      }
//      else                       //NEIGHBOR CELL
//      {
//        const uint64_t cn = face.neighbor_id;
//        const auto& cell_cn = grid.cells[cn];
//        x_cn = cell_cn.centroid;
//
//        const double rho_cn_np1 = U_np1[cn][0];
//
//        sigma_t_cn_np1 = rho_cn_np1 * kappa_t_nph[cn];
//      }
//      const auto x_fcn = x_cn - x_f;
//
//      const double D_cn = -speed_of_light_cmpsh / (3 * sigma_t_cn_np1);
//      const double k_cn = D_cn/x_fcn.Norm();
//
//      const Vec3 x_ccn = x_cn - x_c;
//      const Vec3 kf_np1 = (k_c * k_cn / (k_c + k_cn)) * x_ccn.Normalized();
//
//      const double coeff_LHS = (theta1 / V_c) * A_f.Dot(kf_np1);
//
//      if (not face.has_neighbor) //DEFAULT REFLECTING
//      {
//        //J_f = 0 therefor no connectivity elements
//      }
//      else                       //NEIGHBOR CELL
//      {
//        A[c][c] += -coeff_LHS;
//        const uint64_t cn = face.neighbor_id;
//        A[c][cn] += coeff_LHS;
//      }
//    }//for f in connectivity

    //=========================================== Diagonal and rhs
    A[c][c] += tau + k1 + k4*k5;
    b[c]    += -k3 - k4*k6 + tau*rad_E_c_nphstar - theta2*grad_dot_J_n;
  }//for cell
}