#include "radtranGreyDiff.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

typedef chi_math::VectorN<5> UVector;
typedef std::vector<UVector> GradUTensor;

//###################################################################
/**Computes the gradient of U for the given field data.*/
std::vector<GradUTensor> chi_radhydro::RadTranGreyDiffusion::
ComputeUGradients(const std::vector<UVector>& U) const
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
  // The double minmod
  for (const auto& cell : grid->local_cells)
  {
    const double alpha = 2.0;
    const uint64_t c   = cell.local_id;
    const Vec3&    x_c = cell.centroid;

    const auto& U_c = U[c];

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

      grad_U_c[0] = MinModU({grad_U_c[0], alpha * dU_ds_x});
      grad_U_c[1] = MinModU({grad_U_c[1], alpha * dU_ds_y});
      grad_U_c[2] = MinModU({grad_U_c[2], alpha * dU_ds_z});
    }//for face
  }//for cell

  return grad_U;
}

//###################################################################
/**Computes the gradient of rad_E.*/
std::vector<chi_mesh::Vector3> chi_radhydro::RadTranGreyDiffusion::
ComputeRadEGradients(const std::vector<double> &E) const
{
  std::vector<Vec3> grad_rad_E(num_nodes_local);
  for (const auto& cell_c : grid->local_cells)
  {
    const uint64_t c = cell_c.local_id;
    const auto& fv_view = fv->MapFeView(c);
    const double V_c = fv_view->volume;
    const Vec3&  x_c = cell_c.centroid;

    Vec3 grad_rad_E_c;

    const size_t num_faces = cell_c.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const Vec3   A_f = fv_view->face_area[f]*cell_c.faces[f].normal;
      const Vec3&  x_f = cell_c.faces[f].centroid;

      double weighted_rad_E;

      if (not cell_c.faces[f].has_neighbor) //BOUNDARY CONDITION
      {
        weighted_rad_E = E[c];
      }
      else                                  //NEIGHBOR CELL
      {
        const uint64_t cn = cell_c.faces[f].neighbor_id;
        const Vec3& x_cn = grid->cells[cn].centroid;

        const Vec3 x_c_cn = x_cn - x_c;
        const double d_c_cn = x_c_cn.Norm();

        const double w_c  = (x_f - x_c).Norm()/d_c_cn;
        const double w_cn = (x_cn - x_f).Norm()/d_c_cn;

        weighted_rad_E = w_c*E[c] + w_cn*E[cn];
      }

      grad_rad_E_c = (1/V_c) * A_f * weighted_rad_E;
    }//for face

    grad_rad_E[c] = grad_rad_E_c;
  }//for cell

  //======================================== Apply slope limiting
  // The double minmod
  for (const auto& cell : grid->local_cells)
  {
    const double alpha = 2.0;
    const uint64_t c   = cell.local_id;
    const Vec3&    x_c = cell.centroid;

    const auto& rad_E_c = E[c];

    Vec3& grad_rad_E_c = grad_rad_E[c];

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      double  rad_E_cn;
      Vec3    x_cn;
      if (not face.has_neighbor) //BOUNDARY
      {
        rad_E_cn = E[c];
        x_cn = x_c + 2.0*(face.centroid - x_c);
      }
      else
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid->cells[cn];

        rad_E_cn = E[cn];
        x_cn = cell_cn.centroid;
      }

      const auto dE = rad_E_cn - rad_E_c;
      const auto dx = x_cn - x_c;

      const auto dE_ds_x = (dE/dx.NormSquare())*dx.x;
      const auto dE_ds_y = (dE/dx.NormSquare())*dx.y;
      const auto dE_ds_z = (dE/dx.NormSquare())*dx.z;

      grad_rad_E_c.x = MinMod({grad_rad_E_c[0], alpha * dE_ds_x});
      grad_rad_E_c.y = MinMod({grad_rad_E_c[1], alpha * dE_ds_y});
      grad_rad_E_c.z = MinMod({grad_rad_E_c[2], alpha * dE_ds_z});
    }//for face
  }//for cell

  return grad_rad_E;
}