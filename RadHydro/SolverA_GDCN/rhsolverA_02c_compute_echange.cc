#include "rhsolverA.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

namespace chi_radhydro
{

  SolverA_GDCN::SystemEnergy SolverA_GDCN::
  ComputeSysEnergyChange(
                         const double dt,
                         const std::vector<UVector> &U,
                         const std::vector<double>& rad_E,
                         const std::vector<GradUTensor>& grad_U,
                         const std::vector<chi_mesh::Vector3>& grad_rad_E,
                         const std::map<uint64_t, BCSetting>& bc_setttings)
{
  double Emat = 0.0;
  double Erad = 0.0;
  double mat_adv = 0.0;
  double rad_adv = 0.0;
  //======================================== Determine current amount of energy
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t c         = cell.local_id;
    const auto&    fv_view_c = fv->MapFeView(c);
    const double   V_c       = fv_view_c->volume;

    Emat += V_c * U[c][MAT_E];
    Erad += V_c * rad_E[c];
  }

  //======================================== Determine the net energy entering
  //                                         or leaving the system
  for (const uint64_t c : m_bndry_cells_local_ids)
  {
    const auto& cell       = grid->local_cells[c];
    const auto& fv_view_c  = fv->MapFeView(c);
    const auto& face_areas = fv_view_c->face_area;
    const Vec3& x_cc = cell.centroid;

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const Vec3&  n_f = face.normal;
      const double A_f = face_areas[f];
      const Vec3&  x_fc = face.centroid;

      GradUTensor grad_U_c(3);
      Vec3        grad_rad_E_c;
      if (not grad_U.empty())     grad_U_c     = grad_U[c];
      if (not grad_rad_E.empty()) grad_rad_E_c = grad_rad_E[c];

      const auto U_f     = UplusDXGradU(U[c], x_fc - x_cc, grad_U_c);
      const auto rad_E_f = rad_E[c] + (x_fc - x_cc).Dot(grad_rad_E_c);
      const Vec3   u_f   = chi_radhydro::VelocityFromCellU(U_f);
      const double p_f   = IdealGasPressureFromCellU(U_f, m_gamma);
      const double E_f   = U_f[MAT_E];

      const double u_f_dot_n_f = u_f.Dot(n_f);

      if (not face.has_neighbor) //Boundary face
      {
        if (u_f_dot_n_f > 0.0) //Outgoing
        {
          mat_adv -= dt * A_f * n_f.Dot((E_f + p_f) * u_f);
          rad_adv -= dt * A_f * n_f.Dot((4/3.0) * rad_E[c] * u_f);
        }
        else
        {
          const auto& bc_setting = bc_setttings.at(face.neighbor_id);
          const auto U_upw = MakeUFromBC(bc_setting, U_f);
          const auto rad_E_upw = MakeRadEFromBC(bc_setting, rad_E_f);

          const double E_upw = U_upw[MAT_E];
          const double p_upw = IdealGasPressureFromCellU(U_upw, m_gamma);
          const Vec3   u_upw = VelocityFromCellU(U_upw);

          mat_adv -= dt * A_f * n_f.Dot((E_upw + p_upw)*u_upw);
          rad_adv -= dt * A_f * n_f.Dot((4/3.0)*rad_E_upw*u_upw);
        }
      }
    }//for f
  }//for cell

//  chi::log.Log()
//  << "Emat   : " << Emat << "\n"
//  << "Erad   : " << Erad << "\n"
//  << "mat_adv: " << mat_adv << "\n"
//  << "rad_adv: " << rad_adv << "\n";

  return {Emat, Erad, mat_adv, rad_adv};
}

}//namespace chi_radhydro