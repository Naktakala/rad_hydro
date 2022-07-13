#include "rhsolverA.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

double chi_radhydro::SolverA_GDCN::
ComputeGradDotJ(const chi_mesh::MeshContinuum& grid,
                std::shared_ptr<SDM_FV>&       fv_ref,
                const chi_mesh::Cell&          cell_c,
                const double                   sigma_t_c_n,
                const std::vector<double>&     kappa_t_n,
                const std::vector<UVector>&    U_n,
                const std::vector<double>&     rad_E_n)
{
  const uint64_t c = cell_c.local_id;
  const size_t   num_faces = cell_c.faces.size();

  const auto&    fv_view_c   = fv_ref->MapFeView(c);
  const auto&    face_areas  = fv_view_c->face_area;
  const Vec3&    x_c         = cell_c.centroid;
  const double   V_c         = fv_view_c->volume;

  const double rad_E_c_n = rad_E_n[c];

  double grad_dot_J_n = 0.0;
  for (size_t f=0; f<num_faces; ++f)
  {
    const auto& face = cell_c.faces[f];
    const auto& x_f  = face.centroid;
    const auto  A_f  = face_areas[f] * face.normal;

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

    const double sigma_tf_n   = (sigma_t_c_n + sigma_t_cn_n) / 2.0;

    const double Df_n   = -speed_of_light_cmpsh / (3.0 * sigma_tf_n);

    const Vec3 x_ccn = x_cn - x_c;
    const Vec3 kf_n   = Df_n   * x_ccn / x_ccn.NormSquare();

    const double coeff_RHS = (1.0 / V_c) * A_f.Dot(kf_n);

    if (not face.has_neighbor) //DEFAULT REFLECTING
      grad_dot_J_n += 0.0;
    else                       //NEIGHBOR CELL
      grad_dot_J_n += coeff_RHS * (rad_E_cn_n - rad_E_c_n);

  }//for f in connectivity
//  for (size_t f=0; f<num_faces; ++f)
//  {
//    const auto& face = cell_c.faces[f];
//    const auto& x_f  = face.centroid;
//    const auto  A_f  = face_areas[f] * face.normal;
//    const auto  x_cf = x_f - x_c;
//
//    const double D_c = -speed_of_light_cmpsh / (3 * sigma_t_c_n);
//    const double k_c = D_c/x_cf.Norm();
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
//    const auto x_fcn = x_cn - x_f;
//
//    const double D_cn = -speed_of_light_cmpsh / (3 * sigma_t_cn_n);
//    const double k_cn = D_cn/x_fcn.Norm();
//
//    const Vec3 x_ccn = x_cn - x_c;
//    const Vec3 kf_n  = (k_c * k_cn / (k_c + k_cn)) * x_ccn.Normalized();
//
//    const double coeff_RHS = (1.0 / V_c) * A_f.Dot(kf_n);
//
//    if (not face.has_neighbor) //DEFAULT REFLECTING
//      grad_dot_J_n += 0.0;
//    else                       //NEIGHBOR CELL
//      grad_dot_J_n += coeff_RHS * (rad_E_cn_n - rad_E_c_n);
//
//  }//for f in connectivity

  return grad_dot_J_n;
}