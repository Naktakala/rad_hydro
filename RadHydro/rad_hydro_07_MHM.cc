#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

//###################################################################
/**Generic corrector step of the MHM method.*/
void chi_radhydro::
  MHM_HydroRadEPredictor(const chi_mesh::MeshContinuum& grid,
                         chi_math::SpatialDiscretization_FV& fv,
                         const std::vector<double>&            gamma,
                         double tau,
                         const std::vector<UVector>&           U_a,
                         const std::vector<GradUTensor>&       grad_U_a,
                         const std::vector<double>&            rad_E_a,
                         const std::vector<chi_mesh::Vector3>& grad_rad_E_a,
                         std::vector<UVector>&                 U_a_star,
                         std::vector<double>&                  rad_E_a_star)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector U_c_a    = U_a[c];
    const double  rad_E_c_a = rad_E_a[c];

    const double p_c_n = IdealGasPressureFromCellU(U_c_a, gamma[c]);

    UVector U_c_a_star     = U_c_a;
    double  rad_E_c_a_star = rad_E_c_a;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const double A_f  = fv_view->face_area[f];
      const Vec3&  n_f  = cell.faces[f].normal;
      const Vec3&  x_fc = cell.faces[f].centroid;

      const UVector  U_f   = UplusDXGradU(U_c_a, x_fc - x_cc, grad_U_a[c]);
      const Vec3&    u_f   = Vec3(U_f[1],U_f[2],U_f[3])/U_f[0]; // rho*u/rho

      const UVector F_f = MakeF(U_f, p_c_n, n_f);
      const double rad_E_f = rad_E_c_a + (x_fc - x_cc).Dot(grad_rad_E_a[c]);

      U_c_a_star     -= (1 / tau) * (1 / V_c) * A_f * F_f;
      rad_E_c_a_star -= (1 / tau) * (1 / V_c) * (4.0 / 3) * A_f * n_f.Dot(rad_E_f * u_f);
    }//for f

    U_a_star[c]     = U_c_a_star;
    rad_E_a_star[c] = rad_E_c_a_star;
  }//for cell
}

//###################################################################
/**Generic corrector method for the MHM method.*/
void chi_radhydro::
  MHM_HydroRadECorrector(
  const chi_mesh::MeshContinuum&        grid,
  chi_math::SpatialDiscretization_FV&   fv,
  const std::map<uint64_t, BCSetting>& bc_setttings,
  const std::vector<double>&            gamma,
  double                                tau,
  const std::vector<UVector>&           U_a,
  const std::vector<UVector>&           U_b,
  const std::vector<GradUTensor>&       grad_U_b,
  const std::vector<double>&            rad_E_a,
  const std::vector<double>&            rad_E_b,
  const std::vector<chi_mesh::Vector3>& grad_rad_E_b,
  std::vector<UVector>&                 U_b_star,
  std::vector<double>&                  rad_E_b_star
)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector& U_c_a   = U_a[c];
    const double   rad_E_c_a = rad_E_a[c];

    UVector  U_c_b_star     = U_c_a;
    double   rad_E_c_b_star = rad_E_c_a;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&  face  = cell.faces[f];
      const uint64_t bid = face.neighbor_id;
      const double A_f   = fv_view->face_area[f];
      const Vec3&  n_f   = face.normal;
      const Vec3&  x_fc  = face.centroid;

      const UVector& U_L     = UplusDXGradU(U_b[c], x_fc - x_cc, grad_U_b[c]);
      const double   gamma_L = gamma[c];
      const double   p_L     = IdealGasPressureFromCellU(U_L, gamma_L);

      UVector  U_R     = U_L;
      double   gamma_R = gamma_L;
      double   p_R     = p_L;

      double rad_E_c_f  = rad_E_c_a + (x_fc - x_cc).Dot(grad_rad_E_b[c]);
      Vec3   u_c_f      = VelocityFromCellU(U_L);
      double rad_E_cn_f = rad_E_c_f;
      Vec3   u_cn_f     = u_c_f;

      if (not face.has_neighbor)
      {
        U_R        = MakeUFromBC(bc_setttings.at(bid), U_L);
        rad_E_cn_f = MakeRadEFromBC(bc_setttings.at(bid), rad_E_c_f);
        p_R        = IdealGasPressureFromCellU(U_R, gamma_R);
        u_cn_f     = VelocityFromCellU(U_R);
      }
      else
      {
        const uint64_t cn            = cell.faces[f].neighbor_id;
        const auto&    adjacent_cell = grid.cells[cn];
        const Vec3&    adj_x_cc       = adjacent_cell.centroid;

        U_R     = UplusDXGradU(U_b[cn], x_fc - adj_x_cc, grad_U_b[cn]);
        gamma_R = gamma[cn];
        p_R     = IdealGasPressureFromCellU(U_R, gamma_R);
        rad_E_cn_f = rad_E_b[cn] + (x_fc - adj_x_cc).Dot(grad_rad_E_b[cn]);
        u_cn_f = VelocityFromCellU(U_R);
      }

      //Upwinding rad_Eu
      Vec3 rad_Eu_upw;
      if (u_c_f.Dot(n_f) > 0 and u_cn_f.Dot(n_f) > 0)
        rad_Eu_upw = rad_E_c_f * u_c_f;
      else if (u_c_f.Dot(n_f) < 0 and u_cn_f.Dot(n_f) < 0)
        rad_Eu_upw = rad_E_cn_f * u_cn_f;
      else if (u_c_f.Dot(n_f) > 0 and u_cn_f.Dot(n_f) < 0)
        rad_Eu_upw = rad_E_c_f * u_c_f + rad_E_cn_f * u_cn_f;
      else
        rad_Eu_upw = Vec3(0,0,0);
//      Vec3 rad_Eu_upw = rad_E_c_f * u_c_f;

      const FVector F_hllc_f = HLLC_RiemannSolve(U_L    , U_R,
                                                 p_L    , p_R,
                                                 gamma_L, gamma_R,
                                                 n_f);

      U_c_b_star     -= (1/tau) * (1/V_c) * A_f * F_hllc_f;
      rad_E_c_b_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_Eu_upw);
    }//for f
    U_b_star[c] = U_c_b_star;
    rad_E_b_star[c] = rad_E_c_b_star;
  }//for cell
}

//###################################################################
/**Generic update of density and momentum from lagged
 * radiation momentum deposition.*/
void chi_radhydro::
  DensityMomentumUpdateWithRadMom(
    const chi_mesh::MeshContinuum&      grid,
    chi_math::SpatialDiscretization_FV& fv,
    const std::map<uint64_t, BCSetting>& bc_setttings,
    const std::vector<double>&          kappa_t,
    double tau,
    const std::vector<UVector>&           U_old,
    const std::vector<UVector>&           U_int,
    const std::vector<double>&            rad_E_old,
    std::vector<UVector>&                 U_new
    )
{
  U_new = U_int;
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;

    const double rho_c_old   = U_old[c][0];
    const double rad_E_c_old = rad_E_old[c];

    UVector U_c_new_01 = U_new[c];

    if (std::fabs(kappa_t[c]) < 1.0e-10) continue;

    const double D_c     = -speed_of_light_cmpsh / (3 * rho_c_old * kappa_t[c]);

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&    face = cell.faces[f];
      const uint64_t bid  = face.neighbor_id;
      const Vec3     A_f  = fv_view->face_area[f] * face.normal;
      const Vec3     x_cf = face.centroid - cell.centroid;

      const double k_c = D_c/x_cf.Norm();

      double k_cn = k_c;
      double rad_E_cn_old = rad_E_c_old;
      if (not face.has_neighbor)
        rad_E_cn_old = MakeRadEFromBC(bc_setttings.at(bid), rad_E_c_old);
      else
      {
        const uint64_t cn = face.neighbor_id;
        const auto& adj_cell = grid.cells[cn];
        const Vec3   x_fcn = adj_cell.centroid - cell.faces[f].centroid;

        const double rho_cn_old = U_old[cn][0];
        rad_E_cn_old = rad_E_old[cn];

        const double D_cn = -speed_of_light_cmpsh / (3 * rho_cn_old * kappa_t[cn]);

        k_cn = D_cn/x_fcn.Norm();
      }

      const double rad_E_f = (k_cn * rad_E_cn_old + k_c * rad_E_c_old) /
                             (k_cn + k_c);
      U_c_new_01(1) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.x * rad_E_f;
      U_c_new_01(2) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.y * rad_E_f;
      U_c_new_01(3) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.z * rad_E_f;
    }//for f

    U_new[c] = U_c_new_01;
  }//for cell
}