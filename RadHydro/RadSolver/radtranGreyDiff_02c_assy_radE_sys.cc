#include "radtranGreyDiff.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

void chi_radhydro::RadTranGreyDiffusion::
  AssemblePredictorRadESystem(MatDbl &A, VecDbl &b,
                              chi_mesh::MeshContinuum&  grid_ref,
                              std::shared_ptr<SDM_FV>&  fv_ref,
                              std::vector<UVector>&     U_n,
                              std::vector<UVector>&     U_nstar,
                              std::vector<UVector>&     U_nph,
                              std::vector<double>&      rad_E_n,
                              std::vector<double>&      rad_E_nstar,
                              std::vector<double>&      kappa_t,
                              std::vector<double>&      kappa_a,
                              std::vector<GradUTensor>& grad_U_n,
                              std::vector<Vec3>&        grad_rad_E,
                              std::vector<double>&      Cv,
                              std::vector<double>&      k5_vec,
                              std::vector<double>&      k6_vec,
                              double                    delta_t)
{
  const auto& grid    = grid_ref;
  const auto& fv      = fv_ref;

  k5_vec.assign(grid.local_cells.size(), 0.0);
  k6_vec.assign(grid.local_cells.size(), 0.0);

  for (const auto& cell_c : grid.local_cells)
  {
    const uint64_t c         = cell_c.local_id;
    const size_t   num_faces = cell_c.faces.size();

    const auto&    fv_view_c = fv->MapFeView(c);
    const auto&    face_areas = fv_view_c->face_area;
    const Vec3&    x_c       = cell_c.centroid;
    const double   V_c       = fv_view_c->volume;
    const double   rho_c_nph = U_nph[c][0];
    const double   T_c_n     = IdealGasTemperatureFromCellU(U_n[c],Cv[c]);
    const double   T_c_nstar = IdealGasTemperatureFromCellU(U_nstar[c],Cv[c]);

    const double   E_c_nstar = U_nstar[c][4];

    const double   e_c_nstar   = InternalEnergyFromCellU(U_nstar[c]);
    const double rad_E_c_n     = rad_E_n[c];
    const double rad_E_c_nstar = rad_E_nstar[c];

    const Vec3 u_c_nph = VelocityFromCellU(U_nph[c]);

    const double sigma_tc_nph = rho_c_nph * kappa_t[c];
    const double sigma_ac_nph = rho_c_nph * kappa_a[c];

    const double tau = 1/delta_t;

    double third_grad_rad_E_dot_u = 0.0;
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell_c.faces[f];
      const auto& x_f  = face.centroid;
      const auto  x_fc = x_f - x_c;
      const Vec3 A_f = face_areas[f] * cell_c.faces[f].normal;

      const double  rad_E_f = rad_E_c_n + (x_f - x_c).Dot(grad_rad_E[c]);
      const UVector U_f     = UplusDXGradU(U_n[c], x_fc, grad_U_n[c]);
      const Vec3    u_f     = VelocityFromCellU(U_f);

      third_grad_rad_E_dot_u += (1/V_c)*(1.0/3) * A_f.Dot(rad_E_f * u_f);
    }//for f

    const double k1 = 0.5 * sigma_ac_nph * speed_of_light_cmpsh;
    const double k2 = 4 * pow(T_c_nstar,3)/Cv[c];

    const double k3 = k1*rad_E_c_n
                    - k1*a*pow(T_c_nstar,4)
                    + k1*a*e_c_nstar
                    - k1*a*pow(T_c_n,4)
                    - third_grad_rad_E_dot_u;
    const double k4 = -k1*a*k2;

    const double k5 = k1/(tau*rho_c_nph);
    const double k6 = (k3
                    - tau*0.5*rho_c_nph*pow(u_c_nph.Norm(),2)
                    + tau*E_c_nstar)/
                    (tau*rho_c_nph);

    k5_vec[c] = k5;
    k6_vec[c] = k6;

    //======================= Connectivity and grad_J_n
    double grad_J_n = 0.0;
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell_c.faces[f];
      const auto& x_f  = face.centroid;
      const auto  A_f  = face_areas[f] * face.normal;

      double sigma_tcn_nph = sigma_tc_nph;
      Vec3   x_cn          = x_c + 2*(x_f-x_c);
      double rad_E_cn_n    = rad_E_n[c];

      if (face.has_neighbor)
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid.cells[cn];
        const double kappa_tcn = kappa_t[cn];
        const double rho_cn_nph = U_nph[cn][0];

        sigma_tcn_nph = rho_cn_nph * kappa_tcn;
        x_cn = cell_cn.centroid;
        rad_E_cn_n = rad_E_n[cn];
      }

      const double sigma_tf = (sigma_tc_nph + sigma_tcn_nph)/2.0;
      const double Df = -speed_of_light_cmpsh/(3.0*sigma_tf);

      const Vec3 x_ccn = x_cn - x_c;
      const Vec3 kf = Df * x_ccn / x_ccn.NormSquare();

      const double coeff = (1/(2.0*V_c))*A_f.Dot(kf);

      grad_J_n += coeff*(rad_E_cn_n - rad_E_c_n);
      A[c][c] += -coeff;
      if (face.has_neighbor)
      {
        const uint64_t cn = face.neighbor_id;
        A[c][cn] += coeff;
      }
    }//for f

    //============================ Diagonal and rhs
    A[c][c] += tau + k1 + k4*k5;
    b[c]    += -k3 - k4*k6 + tau*rad_E_c_nstar - grad_J_n;
  }//for cell

}