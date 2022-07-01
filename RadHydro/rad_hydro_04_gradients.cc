#include "rad_hydro.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>

typedef chi_math::VectorN<5> UVector;
typedef std::vector<UVector> GradUTensor;

//###################################################################
/**Min mod operator.*/
double chi_radhydro::MinMod(const std::vector<double> &ais,
                            bool verbose/*=false*/)
{
  if (verbose)
  {
    std::stringstream outp;
    outp << "MinMod: {";
    for (auto a_i : ais) outp << a_i << ",";
    outp << "} ";
    outp << "Signs: {";
    for (auto a_i : ais) outp << std::signbit(a_i) << ",";
    outp << "}";

    chi::log.Log() << outp.str();
  }

  if (ais.empty()) return 0.0;

  const auto same_sign = std::signbit(ais.front());

  for (double a_i : ais)
    if (std::signbit(a_i) != same_sign) return 0.0;

  const bool all_negative = bool(same_sign);

  if (not all_negative) return *std::min_element(ais.begin(), ais.end());
  else                  return *std::max_element(ais.begin(), ais.end());


}

//###################################################################
/**Applies minmod limiter to a UVector*/
chi_radhydro::UVector chi_radhydro::MinModU(const std::vector<UVector> &vec_of_U,
                              bool verbose/*=false*/)
{
  if (vec_of_U.empty())
    throw std::logic_error("chi_hydro::CompInFFlow::MinModU used with no list.");

  UVector minmod_val({0,0,0,0,0});

  const size_t num_elements = vec_of_U.size();
  for (size_t i=0; i<5; ++i)
  {
    std::vector<double> minmod_args(num_elements);

    for (size_t k=0; k<num_elements; ++k)
      minmod_args[k] = vec_of_U[k][static_cast<int>(i)];

    minmod_val(static_cast<int>(i)) = MinMod(minmod_args,verbose);
  }

  return minmod_val;
}

//###################################################################
/**Computes the gradient of U for the given field data.*/
std::vector<GradUTensor> chi_radhydro::
  ComputeUGradients(const std::vector<UVector> &U,
                    const chi_mesh::MeshContinuum &grid,
                    chi_math::SpatialDiscretization_FV &fv,
                    const std::map<uint64_t, BCSetting> &bc_settings)
{
  const size_t num_nodes_local = U.size();
  std::vector<GradUTensor> grad_U(num_nodes_local);
  for (const auto& cell_c : grid.local_cells)
  {
    const uint64_t c = cell_c.local_id;
    const auto& fv_view = fv.MapFeView(c);
    const double V_c = fv_view->volume;
    const Vec3&  x_c = cell_c.centroid;

    GradUTensor grad_U_c(3);

    const size_t num_faces = cell_c.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&  face = cell_c.faces[f];
      const Vec3   A_f = fv_view->face_area[f]*cell_c.faces[f].normal;
      const Vec3&  x_f = cell_c.faces[f].centroid;

      UVector weighted_U;

      if (not face.has_neighbor) //BOUNDARY CONDITION
      {
        const uint64_t bid = face.neighbor_id;
        weighted_U = chi_radhydro::MakeUFromBC(bc_settings.at(bid),U[c]);
      }
      else                       //NEIGHBOR CELL
      {
        const uint64_t cn = face.neighbor_id;
        const Vec3& x_cn = grid.cells[cn].centroid;

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
  for (const auto& cell : grid.local_cells)
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
        const uint64_t bid = face.neighbor_id;
        U_cn = chi_radhydro::MakeUFromBC(bc_settings.at(bid),U[c]);

        x_cn = x_c + 2.0*(face.centroid - x_c);
      }
      else                       //NEIGHBOR CELL
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid.cells[cn];

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
std::vector<chi_mesh::Vector3> chi_radhydro::
ComputeRadEGradients(const std::vector<double> &radE,
                     const chi_mesh::MeshContinuum& grid,
                     chi_math::SpatialDiscretization_FV& fv,
                     const std::map<uint64_t, BCSetting> &bc_settings)
{
  const size_t num_nodes_local = radE.size();
  std::vector<Vec3> grad_rad_E(num_nodes_local);
  for (const auto& cell_c : grid.local_cells)
  {
    const uint64_t c = cell_c.local_id;
    const auto& fv_view = fv.MapFeView(c);
    const double V_c = fv_view->volume;
    const Vec3&  x_c = cell_c.centroid;

    Vec3 grad_rad_E_c;

    const size_t num_faces = cell_c.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto&  face = cell_c.faces[f];
      const Vec3   A_f = fv_view->face_area[f]*cell_c.faces[f].normal;
      const Vec3&  x_f = cell_c.faces[f].centroid;

      double weighted_rad_E;

      if (not face.has_neighbor) //BOUNDARY CONDITION
      {
        const uint64_t bid = face.neighbor_id;
        weighted_rad_E = MakeRadEFromBC(bc_settings.at(bid),radE[c]);
      }
      else                       //NEIGHBOR CELL
      {
        const uint64_t cn = cell_c.faces[f].neighbor_id;
        const Vec3& x_cn = grid.cells[cn].centroid;

        const Vec3 x_c_cn = x_cn - x_c;
        const double d_c_cn = x_c_cn.Norm();

        const double w_c  = (x_f - x_c).Norm()/d_c_cn;
        const double w_cn = (x_cn - x_f).Norm()/d_c_cn;

        weighted_rad_E = w_c * radE[c] + w_cn * radE[cn];
      }

      grad_rad_E_c = (1/V_c) * A_f * weighted_rad_E;
    }//for face

    grad_rad_E[c] = grad_rad_E_c;
  }//for cell

  //======================================== Apply slope limiting
  // The double minmod
  for (const auto& cell : grid.local_cells)
  {
    const double alpha = 2.0;
    const uint64_t c   = cell.local_id;
    const Vec3&    x_c = cell.centroid;

    const auto& rad_E_c = radE[c];

    Vec3& grad_rad_E_c = grad_rad_E[c];

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      double  rad_E_cn;
      Vec3    x_cn;
      if (not face.has_neighbor) //BOUNDARY
      {
        const uint64_t bid = face.neighbor_id;
        rad_E_cn = chi_radhydro::MakeRadEFromBC(bc_settings.at(bid),radE[c]);

        x_cn = x_c + 2.0*(face.centroid - x_c);
      }
      else                       //NEIGHBOR CELL
      {
        const uint64_t cn = face.neighbor_id;
        const auto& cell_cn = grid.cells[cn];

        rad_E_cn = radE[cn];
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