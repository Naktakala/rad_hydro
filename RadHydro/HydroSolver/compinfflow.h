#ifndef RADHYDRO_COMPINFFLOW_H
#define RADHYDRO_COMPINFFLOW_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_hydro
{

enum class CoordinateSystem
{
  NO_GEOMETRY_SET  = 0,
  ONED_SLAB        = 1,
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL   = 3,
  TWOD_CARTESIAN   = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6
};

//Compressible Inviscid Fluid-Flow solver
//Comp         In       F     F ow
//CompInFFlow
class CompInFFlow : public chi_physics::Solver
{
protected:
  typedef std::shared_ptr<SpatialDiscretization_FV> SDMFVPtr;
  typedef std::tuple<chi_mesh::LogicalVolume*, std::string, double> LVFieldSetting;

public:
  double       delta_t_max   = 1.0e-8;
  double       C_cfl         = 0.9;
  CoordinateSystem coordinate_system = CoordinateSystem::NO_GEOMETRY_SET;

  std::vector<LVFieldSetting> logvol_field_settings;

  chi_mesh::MeshContinuumPtr grid;
  chi_math::UnknownManager   uk_man_U;
  SDMFVPtr                   fv;

  uint64_t num_nodes_local = 0;
  uint64_t num_dofs_local = 0;

  std::vector<chi_math::VectorN<4>> U_old, U_new;
  std::vector<double> gamma;
  std::vector<double> rho, u, v, w, p, T, e, Cv;

public:
  explicit CompInFFlow(const std::string& text_name);

  void Initialize() override;
  void Execute() override;
};

}//namespace chi_radhydro

#endif //RADHYDRO_COMPINFFLOW_H
