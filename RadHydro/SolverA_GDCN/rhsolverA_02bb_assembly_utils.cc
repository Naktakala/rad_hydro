#include "rhsolverA.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

//###################################################################
/**Computes Emission-absorption source.
 * Sea = sigma_a c (aT^4 - radE)*/
double chi_radhydro::SolverA_GDCN::
  MakeEmAbsSource(const double sigma_a,
                  const double T,
                  const double radE)
{
  return sigma_a * speed_of_light_cmpsh *
         (a * pow(T,4.0) - radE);
}

//###################################################################
/**Computes GradDotJ.*/
double chi_radhydro::SolverA_GDCN::
ComputeGradDotJ(SimRefs&                       sim_refs,
                const chi_mesh::Cell&          cell_c,
                const double                   sigma_t_c_n,
                const std::vector<double>&     kappa_t_n,
                const std::vector<UVector>&    U_n,
                const std::vector<double>&     rad_E_n)
{
  auto& fv         = sim_refs.fv;
  const auto& grid = sim_refs.grid;

  const uint64_t c = cell_c.local_id;
  const size_t   num_faces = cell_c.faces.size();

  const auto&    fv_view_c   = fv.MapFeView(c);
  const auto&    face_areas  = fv_view_c->face_area;
  const Vec3&    x_c         = cell_c.centroid;
  const double   V_c         = fv_view_c->volume;

  const double rad_E_c_n = rad_E_n[c];

  double grad_dot_J_n = 0.0;
//  for (size_t f=0; f<num_faces; ++f)
//  {
//    const auto& face = cell_c.faces[f];
//    const auto& x_f  = face.centroid;
//    const auto  A_f  = face_areas[f] * face.normal;
//
//    double sigma_t_cn_n   = sigma_t_c_n;
//    Vec3   x_cn          = x_c + 2*(x_f-x_c);
//
//    double rad_E_cn_n = rad_E_n[c];
//    if (not face.has_neighbor) //DEFAULT REFLECTING BC for radE
//    {
//      //Nothing to do for reflecting bc
//    }
//    else                       //NEIGHBOR CELL
//    {
//      const uint64_t cn = face.neighbor_id;
//      const auto& cell_cn = grid.cells[cn];
//      x_cn = cell_cn.centroid;
//
//      const double rho_cn_n   = U_n[cn][0];
//
//      sigma_t_cn_n   = rho_cn_n   * kappa_t_n[cn];
//
//      rad_E_cn_n = rad_E_n[cn];
//    }
//
//    const double sigma_tf_n   = (sigma_t_c_n + sigma_t_cn_n) / 2.0;
//
//    const double Df_n   = -speed_of_light_cmpsh / (3.0 * sigma_tf_n);
//
//    const Vec3 x_ccn = x_cn - x_c;
//    const Vec3 kf_n   = Df_n   * x_ccn / x_ccn.NormSquare();
//
//    const double coeff_RHS = (1.0 / V_c) * A_f.Dot(kf_n);
//
//    if (not face.has_neighbor) //DEFAULT REFLECTING
//      grad_dot_J_n += 0.0;
//    else                       //NEIGHBOR CELL
//      grad_dot_J_n += coeff_RHS * (rad_E_cn_n - rad_E_c_n);
//
//  }//for f in connectivity

  for (size_t f=0; f<num_faces; ++f)
  {
    const auto& face = cell_c.faces[f];
    const auto& x_f  = face.centroid;
    const auto  A_f  = face_areas[f] * face.normal;
    const auto  x_cf = x_f - x_c;

    const double D_c = -speed_of_light_cmpsh / (3 * sigma_t_c_n);
    const double k_c = D_c/x_cf.Norm();

    double sigma_t_cn_n   = sigma_t_c_n;
    Vec3   x_cn          = x_c + 2*(x_f-x_c);

    double rad_E_cn_n = rad_E_n[c];
    if (not face.has_neighbor) //DEFAULT REFLECTING BC for radE
    {
      //Nothing to do for reflecting bc
    }
    else                       //NEIGHBOR CELL
    {
      const uint64_t cn = face.neighbor_id;
      const auto& cell_cn = grid.cells[cn];
      x_cn = cell_cn.centroid;

      const double rho_cn_n   = U_n[cn][0];

      sigma_t_cn_n   = rho_cn_n   * kappa_t_n[cn];

      rad_E_cn_n = rad_E_n[cn];
    }
    const auto x_fcn = x_cn - x_f;

    const double D_cn = -speed_of_light_cmpsh / (3 * sigma_t_cn_n);
    const double k_cn = D_cn/x_fcn.Norm();

    const Vec3 x_ccn = x_cn - x_c;
    const Vec3 kf_n  = (k_c * k_cn / (k_c + k_cn)) * x_ccn.Normalized();

    const double coeff_RHS = (1.0 / V_c) * A_f.Dot(kf_n);

    if (not face.has_neighbor) //DEFAULT REFLECTING
      grad_dot_J_n += 0.0;
    else                       //NEIGHBOR CELL
      grad_dot_J_n += coeff_RHS * (rad_E_cn_n - rad_E_c_n);

  }//for f in connectivity

  return grad_dot_J_n;
}

//###################################################################
/**This does what you think it does.*/
double chi_radhydro::SolverA_GDCN::
  Make3rdGradRadEDot_u(SimRefs&                        sim_refs,
                       const chi_mesh::Cell&           cell,
                       double                          V_c,
                       const std::vector<double>&      face_areas,
                       const std::vector<double>&      rad_E,
                       const std::vector<Vec3>&        grad_rad_E,
                       const std::vector<UVector>&     U,
                       const std::vector<GradUTensor>& grad_U)
{
  const auto&  grid        = sim_refs.grid;
  const double Cv          = sim_refs.Cv;
  const auto&  bc_settings = sim_refs.bc_settings;
  const auto&  kappa_s_function = sim_refs.kappa_s_function;
  const auto&  kappa_a_function = sim_refs.kappa_a_function;

  double third_grad_rad_E_dot_u_nph = 0.0;

  const uint64_t c       = cell.local_id;
  const int      mat_id  = cell.material_id;
  const Vec3&    x_cc    = cell.centroid;

  const UVector& U_c     = U[c];
  const double   rho_c   = U_c[RHO];
  const Vec3     u_c     = VelocityFromCellU(U_c);
  const double   rad_E_c = rad_E[c];
  const double   T_c     = IdealGasTemperatureFromCellU(U_c,Cv);

  const size_t num_faces = cell.faces.size();
  for (size_t f=0; f<num_faces; ++f)
  {
    const auto&  face  = cell.faces[f];
    const uint64_t bid = face.neighbor_id;
    const double A_f   = face_areas[f];
    const Vec3&  n_f   = face.normal;
    const Vec3   x_cf  = face.centroid - cell.centroid;

    Vec3    x_fcn;
    UVector U_cn;
    double  rad_E_cn;

    if (not face.has_neighbor)
    {
      x_fcn    = x_cf;
      U_cn     = MakeUFromBC(bc_settings.at(bid), U_c);
      rad_E_cn = MakeRadEFromBC(bc_settings.at(bid), rad_E_c);
    }
    else
    {
      const uint64_t cn       = cell.faces[f].neighbor_id;
      const auto&    adj_cell = grid.cells[cn];
      x_fcn = adj_cell.centroid - cell.faces[f].centroid;

      U_cn     = U[cn];
      rad_E_cn = rad_E[cn];
    }

    const double rho_cn   = U_cn[RHO];
    const double T_cn     = IdealGasTemperatureFromCellU(U_cn, Cv);

    const double T_f = pow(0.5*(pow(T_c,4.0) + pow(T_cn,4.0)),0.25);

    const double kappa_s_f = ComputeKappaFromLua(T_f, mat_id, kappa_s_function);
    const double kappa_a_f = ComputeKappaFromLua(T_f, mat_id, kappa_a_function);

    const double kappa_t_f = kappa_a_f + kappa_s_f;

    const double D_c  = 1.0 / (rho_c  * kappa_t_f);
    const double D_cn = 1.0 / (rho_cn * kappa_t_f);

    const double k_c  = D_c/x_cf.Norm();
    const double k_cn = D_cn/x_fcn.Norm();

    const double rad_E_f = (k_cn * rad_E_cn + k_c * rad_E_c) /
                           (k_cn + k_c);

    third_grad_rad_E_dot_u_nph += (1/V_c) * (1.0/3) * A_f * n_f.Dot(rad_E_f*u_c);
  }//for f

  return third_grad_rad_E_dot_u_nph;
}


//###################################################################
/**Performs the inverse algebra, contained in the whitepaper, to
 * get the value of internal energy.*/
void chi_radhydro::SolverA_GDCN::
  InverseAlgebraFor_E(SimRefs&                   sim_refs,
                      const std::vector<double>& rad_E,
                      std::vector<UVector>&      U)
{
  for (const auto& cell_c : sim_refs.grid.local_cells)
  {
    const uint64_t c     = cell_c.local_id;
    const double k5_c = sim_refs.e_recon_coeff_k5[c];
    const double k6_c = sim_refs.e_recon_coeff_k6[c];

    double rho_c_np1     = U[c][RHO];
    Vec3   u_c_np1       = VelocityFromCellU(U[c]);
    double u_abs_sqr_np1 = u_c_np1.NormSquare();
    double e_c_np1       = k5_c*rad_E[c] + k6_c;

    double& E_c_np1 = U[c](MAT_E);

    E_c_np1 = rho_c_np1*(0.5 * u_abs_sqr_np1 + e_c_np1);
  }//for cell c
}