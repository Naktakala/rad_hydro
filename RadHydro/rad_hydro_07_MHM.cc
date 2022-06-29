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
                         const std::vector<UVector>&           U_n,
                         const std::vector<GradUTensor>&       grad_U_n,
                         const std::vector<double>&            rad_E_n,
                         const std::vector<chi_mesh::Vector3>& grad_rad_E_n,
                         std::vector<UVector>&                 U_n_star,
                         std::vector<double>&                  rad_E_n_star)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector U_c_n    = U_n[c];
    const double  rad_E_c_n = rad_E_n[c];

    const double p_c_n = IdealGasPressureFromCellU(U_c_n, gamma[c]);

    UVector U_c_n_star     = U_c_n;
    double  rad_E_c_n_star = rad_E_c_n;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const double A_f  = fv_view->face_area[f];
      const Vec3&  n_f  = cell.faces[f].normal;
      const Vec3&  x_fc = cell.faces[f].centroid;

      const UVector  U_f   = UplusDXGradU(U_c_n, x_fc-x_cc, grad_U_n[c]);
      const Vec3&    u_f   = Vec3(U_f[1],U_f[2],U_f[3])/U_f[0]; // rho*u/rho

      const UVector F_f = MakeF(U_f, p_c_n, n_f);
      const double rad_E_f = rad_E_c_n + (x_fc-x_cc).Dot(grad_rad_E_n[c]);

      U_c_n_star     -= (1/tau) * (1/V_c) *           A_f * F_f;
      rad_E_c_n_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_E_f * u_f);
    }//for f

    U_n_star[c]     = U_c_n_star;
    rad_E_n_star[c] = rad_E_c_n_star;
  }//for cell
}

//###################################################################
/**Generic corrector method for the MHM method.*/
void chi_radhydro::
  MHM_HydroRadECorrector(
  const chi_mesh::MeshContinuum&        grid,
  chi_math::SpatialDiscretization_FV&   fv,
  const std::map<int, BCSetting>& bc_setttings,
  const std::vector<double>&            gamma,
  double                                tau,
  const std::vector<UVector>&           U_n,
  const std::vector<UVector>&           U_nph,
  const std::vector<GradUTensor>&       grad_U_nph,
  const std::vector<double>&            rad_E_n,
  const std::vector<double>&            rad_E_nph,
  const std::vector<chi_mesh::Vector3>& grad_rad_E_nph,
  std::vector<UVector>&                 U_nph_star,
  std::vector<double>&                  rad_E_nph_star
)
{
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector& U_c_n   = U_n[c];
    const double   rad_E_c_n = rad_E_n[c];

    UVector  U_c_nph_star     = U_c_n;
    double   rad_E_c_nph_star = rad_E_c_n;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const double A_f  = fv_view->face_area[f];
      const Vec3&  n_f  = cell.faces[f].normal;
      const Vec3&  x_fc = cell.faces[f].centroid;

      const UVector& U_L     = UplusDXGradU(U_nph[c], x_fc-x_cc, grad_U_nph[c]);
      const double   gamma_L = gamma[c];
      const double   p_L     = IdealGasPressureFromCellU(U_L, gamma_L);

      UVector  U_R     = U_L;
      double   gamma_R = gamma_L;
      double   p_R     = p_L;

      double rad_E_c_f  = rad_E_c_n + (x_fc-x_cc).Dot(grad_rad_E_nph[c]);
      Vec3   u_c_f      = chi_radhydro::VelocityFromCellU(U_L);
      double rad_E_cn_f = rad_E_c_f;
      Vec3   u_cn_f     = u_c_f;

      if (not cell.faces[f].has_neighbor)
      {
        const int bid = static_cast<int>(cell.faces[f].neighbor_id);
        if (bc_setttings.count(bid) > 0)
        {
          U_R = MakeUFromBC(bc_setttings.at(bid), U_L);
          p_R = IdealGasPressureFromCellU(U_R, gamma_R);
          rad_E_cn_f = MakeRadEFromBC(bc_setttings.at(bid), rad_E_c_f);
          u_cn_f = chi_radhydro::VelocityFromCellU(U_R);
        }
      }
      else
      {
        const uint64_t cn            = cell.faces[f].neighbor_id;
        const auto&    adjacent_cell = grid.cells[cn];
        const Vec3&    adj_x_cc       = adjacent_cell.centroid;

        U_R     = UplusDXGradU(U_nph[cn], x_fc-adj_x_cc, grad_U_nph[cn]);
        gamma_R = gamma[cn];
        p_R     = IdealGasPressureFromCellU(U_R, gamma_R);
        rad_E_cn_f = rad_E_nph[cn] + (x_fc-adj_x_cc).Dot(grad_rad_E_nph[cn]);
        u_cn_f = chi_radhydro::VelocityFromCellU(U_R);
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

      const FVector F_hllc_f = HLLC_RiemannSolve(U_L    , U_R,
                                                 p_L    , p_R,
                                                 gamma_L, gamma_R,
                                                 n_f);

      U_c_nph_star -= (1/tau) * (1/V_c)* A_f * F_hllc_f;
      rad_E_c_nph_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_Eu_upw);
    }//for f
    U_nph_star[c] = U_c_nph_star;
    rad_E_nph_star[c] = rad_E_c_nph_star;
  }//for cell
}

//###################################################################
/**Generic update of density and momentum from lagged
 * radiation momentum deposition.*/
void chi_radhydro::
  DensityMomentumUpdateWithRadMom(
    const chi_mesh::MeshContinuum&      grid,
    chi_math::SpatialDiscretization_FV& fv,
    const std::map<int, BCSetting>& bc_setttings,
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
      const Vec3 A_f  = fv_view->face_area[f] * cell.faces[f].normal;
      const Vec3 x_cf = cell.faces[f].centroid - cell.centroid;

      const double k_c = D_c/x_cf.Norm();

      double k_cn = k_c;
      double rad_E_cn_old = rad_E_c_old;
      if (not cell.faces[f].has_neighbor)
      {
        const int bid = static_cast<int>(cell.faces[f].neighbor_id);
        if (bc_setttings.count(bid) > 0)
          rad_E_cn_old = MakeRadEFromBC(bc_setttings.at(bid), rad_E_c_old);
        else
          rad_E_cn_old = rad_E_c_old;
      }
      if (cell.faces[f].has_neighbor)
      {
        const uint64_t cn = cell.faces[f].neighbor_id;
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