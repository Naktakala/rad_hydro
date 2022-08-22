#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

//###################################################################
/**Generic corrector step of the MHM method.*/
void chi_radhydro::
  MHM_HydroPredictor(SimRefs&                              sim_refs,
                     const double                          tau,
                     const std::vector<UVector>&           U_a,
                     const std::vector<GradUTensor>&       grad_U_a,
                     std::vector<UVector>&                 U_a_star)
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

    const double p_c_n = IdealGasPressureFromCellU(U_c_a, gamma);

    UVector U_c_a_star     = U_c_a;

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

      U_c_a_star     -= (1/tau)*(1/V_c)*A_f * F_f;
    }//for f

    U_a_star[c]     = U_c_a_star;
  }//for cell
}

//###################################################################
/**Generic corrector method for the MHM method.*/
void chi_radhydro::
  MHM_HydroCorrector(
  SimRefs&                              sim_refs,
  double                                tau,
  const std::vector<UVector>&           U_old,
  const std::vector<UVector>&           U_int,
  const std::vector<GradUTensor>&       grad_U_int,
  std::vector<UVector>&                 U_int_star
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

    const auto& grad_U_c_int     = grad_U_int[c];

    UVector  U_c_b_star     = U_old[c];

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&  face  = cell.faces[f];
      const uint64_t bid = face.neighbor_id;
      const double A_f   = fv_view->face_area[f];
      const Vec3&  n_f   = face.normal;
      const Vec3&  x_fc  = face.centroid;

      const UVector& U_L       = UplusDXGradU(U_c_int, x_fc - x_cc, grad_U_c_int);

      UVector U_R;

      if (not face.has_neighbor)
      {
        U_R        = MakeUFromBC(bc_settings.at(bid), U_c_int);
      }
      else
      {
        const uint64_t cn            = cell.faces[f].neighbor_id;
        const auto&    adjacent_cell = grid.cells[cn];
        const Vec3&    x_cn          = adjacent_cell.centroid;

        U_R        = UplusDXGradU(U_int[cn], x_fc - x_cn, grad_U_int[cn]);
      }

      const FVector F_hllc_f = HLLC_RiemannSolve(U_L, U_R, gamma, n_f);

      U_c_b_star     -= (1/tau) * (1/V_c) * A_f * F_hllc_f;
    }//for f

    U_int_star[c] = U_c_b_star;
  }//for cell
}