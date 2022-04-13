#include "compinfflow.h"

#include "chi_log.h"
extern ChiLog& chi_log;

typedef chi_math::VectorN<5> UVector;
typedef std::vector<UVector> GradUTensor;

//###################################################################
/**Computes the gradient of U for the given field data.*/
std::vector<GradUTensor> chi_hydro::CompInFFlow::
  ComputeGradients(const std::vector<UVector>& U) const
{
  std::vector<GradUTensor> grad_U(num_nodes_local);
  for (const auto& cell_c : grid->local_cells)
  {
    const uint64_t c = cell_c.local_id;
    const auto& fv_view = fv->MapFeView(c);
    const double V_c = fv_view->volume;
    const Vec3&  x_c = cell_c.centroid;

    GradUTensor grad_U_c(3);

    const size_t num_faces = cell_c.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const Vec3   A_f = fv_view->face_area[f]*cell_c.faces[f].normal;
      const Vec3&  x_f = cell_c.faces[f].centroid;

      UVector weighted_U;

      if (not cell_c.faces[f].has_neighbor) //BOUNDARY CONDITION
      {
        weighted_U = U[c];
      }
      else                                  //NEIGHBOR CELL
      {
        const uint64_t cn = cell_c.faces[f].neighbor_id;
        const Vec3& x_cn = grid->cells[cn].centroid;

        const Vec3 x_c_cn = x_cn - x_c;
        const double d_c_cn = x_c_cn.Norm();

        const double w_c  = (x_f - x_c).Norm()/d_c_cn;
        const double w_cn = (x_cn - x_f).Norm()/d_c_cn;

        weighted_U = w_c*U[c] + w_cn*U[cn];
      }

      grad_U_c[0] += (1/V_c) * A_f.x * weighted_U; //dU_dx_c
      grad_U_c[1] += (1/V_c) * A_f.y * weighted_U; //dU_dy_c
      grad_U_c[2] += (1/V_c) * A_f.z * weighted_U; //dU_dz_c
    }//for face

    grad_U[c] = grad_U_c;
  }//for cell

  //======================================== Apply slope limiting
  for (const auto& cell : grid->local_cells) //TODO: Generalize to 2D/3D
  {
    if (cell.Type() != chi_mesh::CellType::SLAB)
      throw std::logic_error("CompInFFlow::ComputeGradients does not support"
                             " cells other than SLAB.");
    const uint64_t c   = cell.local_id;
    const Vec3&    x_c = cell.centroid;

    Vec3 x_cm1;
    Vec3 x_cp1;

    const auto U_c      = U[c];
    UVector    U_cm1;
    UVector    U_cp1;

    GradUTensor& grad_U_c = grad_U[c];

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      UVector U_cn;
      Vec3    x_cn;
      if (not face.has_neighbor) //BOUNDARY
      {
        U_cn = U[c];
        x_cn = x_c + 2.0*(face.centroid - x_c);
      }
      else
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid->cells[cn];

        U_cn = U[cn];
        x_cn = cell_cn.centroid;
      }

      const auto dU = U_cn - U_c;
      const auto dx = x_cn - x_c;

      const auto dU_ds_x = (dU/dx.NormSquare())*dx.x;
      const auto dU_ds_y = (dU/dx.NormSquare())*dx.y;
      const auto dU_ds_z = (dU/dx.NormSquare())*dx.z;

      if (f==0) x_cm1 = x_cn;
      if (f==1) x_cp1 = x_cn;

      if (f==0) U_cm1 = U_cn;
      if (f==1) U_cp1 = U_cn;
    }//for face

    const double dx_cmh = (x_c - x_cm1).Norm();
    const double dx_cph = (x_c - x_cp1).Norm();

    grad_U_c[0] = UVector({0,0,0,0,0});
    grad_U_c[1] = UVector({0,0,0,0,0});
    grad_U_c[2] = MinModU({(U_cm1-U_c)/dx_cmh,
                           (U_cp1-U_cm1)/(dx_cmh+dx_cph),
                           (U_cp1-U_c)/dx_cph});

//    chi_log.Log() << "c: " << c
//                  << " dU_dx: " << PrintU(grad_U_c[0])
//                  << " dU_dy: " << PrintU(grad_U_c[1])
//                  << " dU_dz: " << PrintU(grad_U_c[2]);
  }//for cell



  return grad_U;
}