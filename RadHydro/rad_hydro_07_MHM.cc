#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

//###################################################################
/**Generic corrector step of the MHM method.*/
void chi_radhydro::
  MHM_HydroRadEPredictor(SimRefs&                              sim_refs,
                         const double                          tau,
                         const std::vector<UVector>&           U_a,
                         const std::vector<GradUTensor>&       grad_U_a,
                         const std::vector<double>&            rad_E_a,
                         const std::vector<chi_mesh::Vector3>& grad_rad_E_a,
                         std::vector<UVector>&                 U_a_star,
                         std::vector<double>&                  rad_E_a_star)
{
  const auto&  grid  = sim_refs.grid;
        auto&  fv    = sim_refs.fv;
  const double gamma = sim_refs.gamma;

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector U_c_a    = U_a[c];
    const double  rad_E_c_a = rad_E_a[c];

    const double p_c_n = IdealGasPressureFromCellU(U_c_a, gamma);

    UVector U_c_a_star     = U_c_a;
    double  rad_E_c_a_star = rad_E_c_a;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const double A_f  = fv_view->face_area[f];
      const Vec3&  n_f  = cell.faces[f].normal;
      const Vec3&  x_fc = cell.faces[f].centroid;
      const Vec3   x_fc_cc = x_fc - x_cc;

      const UVector  U_f   = UplusDXGradU(U_c_a, x_fc_cc, grad_U_a[c]);
      const Vec3&    u_f   = VelocityFromCellU(U_f);

      const UVector F_f = MakeF(U_f, p_c_n, n_f);
      const double rad_E_f = rad_E_c_a + x_fc_cc.Dot(grad_rad_E_a[c]);

      U_c_a_star     -= (1/tau)*(1/V_c)*A_f * F_f;
      rad_E_c_a_star -= (1/tau)*(1/V_c)*A_f * (4.0 / 3) * n_f.Dot(rad_E_f * u_f);
    }//for f

    U_a_star[c]     = U_c_a_star;
    rad_E_a_star[c] = rad_E_c_a_star;
  }//for cell
}

//###################################################################
/**Generic corrector method for the MHM method.*/
void chi_radhydro::
  MHM_HydroRadECorrector(
  SimRefs&                              sim_refs,
  double                                tau,
  const std::vector<UVector>&           U_old,
  const std::vector<UVector>&           U_int,
  const std::vector<GradUTensor>&       grad_U_int,
  const std::vector<double>&            rad_E_old,
  const std::vector<double>&            rad_E_int,
  const std::vector<chi_mesh::Vector3>& grad_rad_E_int,
  std::vector<UVector>&                 U_int_star,
  std::vector<double>&                  rad_E_int_star
)
{
  const auto&  grid        = sim_refs.grid;
        auto&  fv          = sim_refs.fv;
  const auto&  bc_settings = sim_refs.bc_settings;
  const double gamma       = sim_refs.gamma;

  for (const auto& cell : grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector& U_c_int       = U_int[c];
    const double&  rad_E_c_int   = rad_E_int[c];

    const auto& grad_U_c_int     = grad_U_int[c];
    const auto& grad_rad_E_c_int = grad_rad_E_int[c];

    UVector  U_c_b_star     = U_old[c];
    double   rad_E_c_b_star = rad_E_old[c];

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&  face  = cell.faces[f];
      const uint64_t bid = face.neighbor_id;
      const double A_f   = fv_view->face_area[f];
      const Vec3&  n_f   = face.normal;
      const Vec3&  x_fc  = face.centroid;

      const UVector& U_L       = UplusDXGradU(U_c_int, x_fc - x_cc, grad_U_c_int);
      const double   rad_E_c_f = rad_E_c_int + (x_fc - x_cc).Dot(grad_rad_E_int[c]);
      const Vec3     u_c_f     = VelocityFromCellU(U_L);

      UVector U_R;
      double  rad_E_cn_f;
      Vec3    u_cn_f;

      if (not face.has_neighbor)
      {
        U_R        = MakeUFromBC(bc_settings.at(bid), U_c_int);
        rad_E_cn_f = MakeRadEFromBC(bc_settings.at(bid), rad_E_c_f);
        u_cn_f     = VelocityFromCellU(U_R);
      }
      else
      {
        const uint64_t cn            = cell.faces[f].neighbor_id;
        const auto&    adjacent_cell = grid.cells[cn];
        const Vec3&    x_cn          = adjacent_cell.centroid;

        U_R        = UplusDXGradU(U_int[cn], x_fc - x_cn, grad_U_int[cn]);
        rad_E_cn_f = rad_E_int[cn] + (x_fc - x_cn).Dot(grad_rad_E_int[cn]);
        u_cn_f     = VelocityFromCellU(U_R);
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

      const FVector F_hllc_f = HLLC_RiemannSolve(U_L, U_R, gamma, n_f);

      U_c_b_star     -= (1/tau) * (1/V_c) * A_f * F_hllc_f;
      rad_E_c_b_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_Eu_upw);
    }//for f

    U_int_star[c] = U_c_b_star;
    rad_E_int_star[c] = rad_E_c_b_star;
  }//for cell
}

//###################################################################
/**Generic update of density and momentum from lagged
 * radiation momentum deposition.*/
void chi_radhydro::
  DensityMomentumUpdateWithRadMom(
    SimRefs&                              sim_refs,
    const std::vector<double>&            kappa_t,
    double                                tau,
    const std::vector<UVector>&           U_old,
    const std::vector<UVector>&           U_int,
    const std::vector<double>&            rad_E_old,
    std::vector<UVector>&                 U_new
    )
{
  U_new = U_int;

  const auto& Cv          = sim_refs.Cv;
  const auto& bc_settings = sim_refs.bc_settings;
  const auto& grid        = sim_refs.grid;
        auto& fv          = sim_refs.fv;

  const auto& kappa_s_function = sim_refs.kappa_s_function;
  const auto& kappa_a_function = sim_refs.kappa_a_function;

  for (const auto& cell : sim_refs.grid.local_cells)
  {
    const uint64_t c       = cell.local_id;
    const int      mat_id  = cell.material_id;
    const auto&    fv_view = fv.MapFeView(c);
    const double   V_c     = fv_view->volume;

    const UVector& U_c_old     = U_old[c];
    const double   rho_c_old   = U_old[c][RHO];
    const double   rad_E_c_old = rad_E_old[c];
    const double   T_c_old     = IdealGasTemperatureFromCellU(U_c_old,Cv);

    UVector U_c_new_01 = U_new[c];

    if (std::fabs(kappa_t[c]) < 1.0e-10) continue;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&    face = cell.faces[f];
      const uint64_t bid  = face.neighbor_id;
      const Vec3     A_f  = fv_view->face_area[f] * face.normal;
      const Vec3     x_cf = face.centroid - cell.centroid;

      Vec3    x_fcn = x_cf;
      UVector U_cn_old;
      double  rad_E_cn_old;

      if (not face.has_neighbor)
      {
        U_cn_old     = MakeUFromBC(bc_settings.at(bid), U_c_old);
        rad_E_cn_old = MakeRadEFromBC(bc_settings.at(bid), rad_E_c_old);
      }
      else
      {
        const uint64_t cn = face.neighbor_id;
        const auto& adj_cell = grid.cells[cn];
        x_fcn = adj_cell.centroid - cell.faces[f].centroid;

        U_cn_old     = U_old[cn];
        rad_E_cn_old = rad_E_old[cn];
      }

      using cdouble = const double;

      cdouble rho_cn_old   = U_cn_old[RHO];
      cdouble T_cn_old     = IdealGasTemperatureFromCellU(U_cn_old, Cv);

      cdouble T_f = pow(0.5*(pow(T_c_old,4.0) + pow(T_cn_old,4.0)),0.25);

      cdouble kappa_s_f = ComputeKappaFromLua(T_f, mat_id, kappa_s_function);
      cdouble kappa_a_f = ComputeKappaFromLua(T_f, mat_id, kappa_a_function);

      cdouble kappa_t_f = kappa_a_f + kappa_s_f;

      cdouble D_c  = - speed_of_light_cmpsh / (3.0*(rho_c_old  * kappa_t_f));
      cdouble D_cn = - speed_of_light_cmpsh / (3.0*(rho_cn_old * kappa_t_f));

      cdouble k_c  = D_c/x_cf.Norm();
      cdouble k_cn = D_cn/x_fcn.Norm();

      cdouble rad_E_f = (k_cn * rad_E_cn_old + k_c * rad_E_c_old) /
                        (k_cn + k_c);

      U_c_new_01(1) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.x * rad_E_f;
      U_c_new_01(2) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.y * rad_E_f;
      U_c_new_01(3) -= (1 / tau) * (1 / V_c) * (1.0 / 3) * A_f.z * rad_E_f;
    }//for f

    U_new[c] = U_c_new_01;
  }//for cell
}