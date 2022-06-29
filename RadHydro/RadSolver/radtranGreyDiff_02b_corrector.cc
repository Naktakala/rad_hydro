#include "radtranGreyDiff.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/chi_math_banded_solvers.h"
#include "radtran_structs.h"

void chi_radhydro::RadTranGreyDiffusion::
  Corrector(double dt,
            std::vector<UVector>&           U_n,
            std::vector<UVector>&           U_nph,
            std::vector<UVector>&           U_np1,
            std::vector<GradUTensor>&       grad_U,

            std::vector<double>&            rad_E_n,
            std::vector<double>&            rad_E_nph,
            std::vector<double>&            rad_E_np1,
            std::vector<chi_mesh::Vector3>& grad_rad_E)
{
  chi::log.Log0() << "Executing corrector";

  std::vector<UVector> U_nph_star(num_nodes_local);
  std::vector<double>  rad_E_nph_star(num_nodes_local);

  //##################################################### CORRECTOR
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
    const double T = chi_radhydro::IdealGasTemperatureFromCellU(U_nph[c], Cv[c]);
    kappa_a[c] = kappa_1/(kappa_2*pow(T,n_exponent) + kappa_3);

    kappa_t[c] = kappa_s + kappa_a[c];
  }//for cell

  //=================================== Advection of U and rad_E
  // Advance full timestep
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv->MapFeView(c);
    const double   V_c     = fv_view->volume;
    const Vec3&    x_cc    = cell.centroid;

    const UVector& U_c_n   = U_n[c];
    const double   rad_E_c_n = rad_E_n[c];

    UVector  U_c_nph_star     = U_c_n;
    double   rad_E_c_nph_star = rad_E_c_n;

    const double   tau     = dt;

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

      double rad_E_c_f  = rad_E_c_n + (x_fc-x_cc).Dot(grad_rad_E[c]);
      Vec3   u_c_f      = chi_radhydro::VelocityFromCellU(U_L);
      double rad_E_cn_f = rad_E_c_f;
      Vec3   u_cn_f     = u_c_f;

      if (cell.faces[f].has_neighbor)
      {
        const uint64_t cn            = cell.faces[f].neighbor_id;
        const auto&    adjacent_cell = grid->cells[cn];
        const Vec3&    adj_x_cc       = adjacent_cell.centroid;

        U_R     = UplusDXGradU(U_nph[cn], x_fc-adj_x_cc, grad_U[cn]);
        gamma_R = gamma[cn];
        p_R     = IdealGasPressureFromCellU(U_R, gamma_R);
        rad_E_cn_f = rad_E_nph[cn] + (x_fc-adj_x_cc).Dot(grad_rad_E[cn]);
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

      U_c_nph_star -= (dt/V_c)* A_f * F_hllc_f;
      rad_E_c_nph_star -= (1/tau) * (1/V_c) * (4.0/3) * A_f * n_f.Dot(rad_Eu_upw);
    }//for f
    U_nph_star[c] = U_c_nph_star;
    rad_E_nph_star[c] = rad_E_c_nph_star;
  }//for cell

  // Update density and momentum to np1
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t c       = cell.local_id;
    const auto&    fv_view = fv->MapFeView(c);
    const double   V_c     = fv_view->volume;

    const double   tau     = 0.5*dt;

    const double D_c     = -speed_of_light_cmpsh / (3 * rho[c] * kappa_t[c]);
    const double rad_E_c = rad_E_nph[c];

    UVector U_c_np1_01 = U_nph_star[c];

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
        rad_E_cn = rad_E_nph[cn];
      }

      const double rad_E_f = (k_cn*rad_E_cn + k_c*rad_E_c)/
                             (k_cn + k_c);
      U_c_np1_01(1) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.x * rad_E_f;
      U_c_np1_01(2) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.y * rad_E_f;
      U_c_np1_01(3) -= (1/tau) * (1/V_c) * (1.0/3) * A_f.z * rad_E_f;
    }//for f

    U_np1[c] = U_c_np1_01;
  }//for cell

  // Internal energy and radiation energy
  {
    std::vector<double> k5,k6;
//    chi_radhydro::DataBlockA data(
//      *grid,
//      fv,
//      U_n,
//      U_n_star,
//      U_nph,
//      rad_E_n,
//      rad_E_n_star,
//      kappa_t,
//      kappa_a,
//      grad_U,
//      grad_rad_E,
//      Cv,
//      k5,k6,
//      dt
//    );
    MatDbl A(num_nodes_local, VecDbl(num_nodes_local,0.0));
    VecDbl b(num_nodes_local, 0.0);

//    AssemblePredictorRadESystem(A, b, data);
    rad_E_nph = chi_math::TDMA(A,b);

    for (const auto& cell_c : grid->local_cells)
    {
      const uint64_t c = cell_c.local_id;
      double rho_c_nph = U_nph[c][0];
      Vec3   u_c_nph   = Vec3(U_nph[c][1],U_nph[c][2],U_nph[c][3])/rho_c_nph;
      double u_abs_sqr = u_c_nph.NormSquare();
      double e_c_nph = k5[c]*rad_E_nph[c] + k6[c];

      double& E_c_nph = U_nph[c](4);

      E_c_nph = 0.5*rho_c_nph*(u_abs_sqr + e_c_nph);
    }//for cell c
  }//scope internal e and rad_E

  chi::log.Log0() << "Done executing corrector";
}