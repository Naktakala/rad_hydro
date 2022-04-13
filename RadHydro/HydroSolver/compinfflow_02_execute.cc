#include "compinfflow.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Executes the solver.*/
void chi_hydro::CompInFFlow::Execute()
{
  chi_log.Log() << "\nExecuting CompInFFlow solver\n\n";

  std::vector<UVector>              U_n(num_nodes_local);
  std::vector<std::vector<UVector>> U_nph(num_nodes_local);
  std::vector<UVector>              U_np1(num_nodes_local);

  FieldsToU(U_n);

  for (int n=0; n<num_timesteps; ++n)
  {
    //================================= Compute delta_t
    const double delta_t_hydro = ComputeCourantLimitDelta_t();
    const double dt            = std::min(delta_t_max, delta_t_hydro);

    chi_log.Log() << "Timestep " << n << " with dt=" << dt;

    //================================= Compute gradients
    const auto grad_U = ComputeGradients(U_n); //Populates grad_U

    //================================= Advance half timestep using gradients
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c       = cell.local_id;
      const auto&    fv_view = fv->MapFeView(c);
      const double   V_c     = fv_view->volume;
      const Vec3&    x_c     = cell.centroid;

      const UVector U_c_n   = U_n[c];

      //========================== Reconstruct edge values
      const size_t num_faces = cell.faces.size();
      std::vector<UVector> U_f(num_faces);
      for (size_t f=0; f<num_faces; ++f)
      {
        const Vec3&  x_f  = cell.faces[f].centroid;
        const Vec3   x_fc = x_f - x_c;

        U_f[f] = U_c_n + x_fc.x * grad_U[c][0] +
                         x_fc.y * grad_U[c][1] +
                         x_fc.z * grad_U[c][2];
      }

      //========================== Advance edge values 0.5 dt
      std::vector<UVector> bar_U_f = U_f;
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f  = fv_view->face_area[f];
        const Vec3&  n_f  = cell.faces[f].normal;

        auto T    = MakeTransformationMatrix(n_f, coordinate_system);
        auto Tinv = T.Inverse();

        UVector F_f = MakeF(T*U_f[f], p[c]);

        for (size_t fp=0; fp<num_faces; ++fp) //fp=fprime
          bar_U_f[fp] -= (0.5*dt / V_c) * A_f * Tinv * F_f;
      }//for f

      U_nph[c] = bar_U_f;
    }//for cell

    //====================================== Execute Riemann solves
    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c = cell.local_id;
      const auto& fv_view = fv->MapFeView(c);
      const double V_c = fv_view->volume;

      const UVector&              U_c_n   = U_n[c];
      const std::vector<UVector>& U_c_nph = U_nph[c];


      UVector U_c_np1 = U_c_n;

      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const double A_f = fv_view->face_area[f];
        const Vec3&  n_f = cell.faces[f].normal;

        const UVector& U_L     = U_c_nph[f];
        const double   gamma_L = gamma[c];
        const double   p_L     = PressureFromCellU(U_L, gamma_L);

        UVector U_R = U_L;
        double gamma_R = gamma_L;
        double p_R = p_L;

        if (cell.faces[f].has_neighbor)
        {
          const uint64_t cn = cell.faces[f].neighbor_id;
          const size_t face_index_mapped = cell_face_mappings[c][f];

          U_R     = U_nph[cn][face_index_mapped];
          gamma_R = gamma[cn];
          p_R     = PressureFromCellU(U_R, gamma_R);
        }

        auto T    = MakeTransformationMatrix(n_f, coordinate_system);
        auto Tinv = T.Inverse();

        bool verbose = false;
        const FVector F_hllc_f = HLLC_RiemannSolve(T*U_L,T*U_R,
                                                   p_L    ,p_R,
                                                   gamma_L,gamma_R,
                                                   verbose);

        U_c_np1 -= (dt/V_c)* A_f * Tinv * F_hllc_f;
      }//for f
      U_np1[c] = U_c_np1;
    }//for cell

    U_n = U_np1;
    UToFields(U_n);

    auto PrintField = [](const std::vector<double>& field,
                         const std::string& field_name)
    {
      std::stringstream out; out << field_name << ": ";
      for (double val : field)
      {
        if (std::fabs(val) > 1.0e-10)
          out << std::setprecision(2) << std::scientific << val << " ";
        else
          out << std::setprecision(2) << std::scientific << 0.0 << " ";
      }

      return out.str();
    };

    std::vector<double> xc(num_nodes_local,0.0);
    for (const auto& cell : grid->local_cells)
      xc[cell.local_id] = cell.centroid.z;

    std::ofstream ofile;
    ofile.open("Output.txt",std::ofstream::out|std::ofstream::app);

    ofile << "Timestep: " <<std::setw(4) << n << " dt=" << dt << " ";
    ofile << PrintField(w  ,"u  ") << "\n";

    ofile.close();
  }//for n timesteps




  chi_log.Log() << "\nDone executing CompInFFlow solver\n\n";
}