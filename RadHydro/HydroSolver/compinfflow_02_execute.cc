#include "compinfflow.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Executes the solver.*/
void chi_hydro::CompInFFlow::Execute()
{
  chi_log.Log() << "\nExecuting CompInFFlow solver\n\n";

  std::vector<UVector>&             U_n = U_old;
  std::vector<UVector>              U_nph(num_nodes_local);
  std::vector<UVector>              U_np1(num_nodes_local);

  FieldsToU(U_n);

  double time = 0.0;
  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double delta_t_hydro = ComputeCourantLimitDelta_t();
    const double dt            = std::min(delta_t_max, delta_t_hydro);

    chi_log.Log() << "Timestep " << n << " with dt=" << dt << " t=" << time;

    //================================= Compute gradients
    //Populates grad_U with slope limiters
    const auto grad_U = ComputeGradients(U_n);

    //================================= Advance half timestep using gradients
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;
      const Vec3&    x_cc    = cell.centroid;

      const UVector U_c_n   = U_n[c];
            UVector U_c_nph = U_c_n;

      //========================== Advance 0.5 dt
      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f  = fv_view->face_area[f];
        const Vec3&  n_f  = cell.faces[f].normal;
        const Vec3&  x_fc = cell.faces[f].centroid;

        const UVector& U_f = UplusDXGradU(U_c_n, x_fc-x_cc, grad_U[c]);

        const UVector F_f = MakeF(U_f, p[c], n_f);

        U_c_nph -= (0.5*dt/V_c) * A_f * F_f;
      }//for f

      U_nph[c] = U_c_nph;
    }//for cell

    //====================================== Execute Riemann solves
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;
      const Vec3&    x_cc    = cell.centroid;

      const UVector& U_c_n   = U_n[c];
            UVector  U_c_np1 = U_c_n;

      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f  = fv_view->face_area[f];
        const Vec3&  n_f  = cell.faces[f].normal;
        const Vec3&  x_fc = cell.faces[f].centroid;

        const UVector& U_L     = UplusDXGradU(U_nph[c], x_fc-x_cc, grad_U[c]);
        const double   gamma_L = gamma[c];
        const double   p_L     = PressureFromCellU(U_L, gamma_L);

              UVector  U_R     = U_L;
              double   gamma_R = gamma_L;
              double   p_R     = p_L;

        if (cell.faces[f].has_neighbor)
        {
          const uint64_t cn            = cell.faces[f].neighbor_id;
          const auto&    adjacent_cell = grid->cells[cn];
          const Vec3&    adj_x_cc       = adjacent_cell.centroid;

          U_R     = UplusDXGradU(U_nph[cn], x_fc-adj_x_cc, grad_U[cn]);
          gamma_R = gamma[cn];
          p_R     = PressureFromCellU(U_R, gamma_R);
        }

        const FVector F_hllc_f = HLLC_RiemannSolve(U_L    , U_R,
                                                   p_L    , p_R,
                                                   gamma_L, gamma_R,
                                                   n_f);

        U_c_np1 -= (dt/V_c)* A_f * F_hllc_f;
      }//for f
      U_np1[c] = U_c_np1;
    }//for cell

    U_n = U_np1;
    UToFields(U_n);

    time += dt;
    if (t_max >= 0.0 and time >= t_max)
      break;
  }//for n timesteps

  chi_log.Log() << "\nDone executing CompInFFlow solver\n\n";
}