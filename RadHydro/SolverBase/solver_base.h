#ifndef RADHYDRO_SOLVER_BASE_H
#define RADHYDRO_SOLVER_BASE_H

#include "rad_hydro.h"
#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMesh/chi_mesh.h"

#include "ChiMath/chi_math.h"

#include <map>

namespace chi_radhydro
{

class RadHydroSolver : public chi_physics::Solver
{
protected:
  typedef std::shared_ptr<chi_mesh::LogicalVolume> LogicalVolumePtr;
  typedef std::tuple<LogicalVolumePtr, std::string, double> LVFieldSetting;
  typedef std::vector<double> ScalarField;
  typedef chi_math::SpatialDiscretization_FV SDM_FV;

protected:
  std::shared_ptr<chi_mesh::MeshContinuum> grid;
  chi_math::UnknownManager ONE_DOF_PER_NODE;
  std::shared_ptr<SDM_FV>                  fv;

public:
  std::vector<LVFieldSetting>        logvol_field_settings;
  std::map<std::string, ScalarField> scalar_fields;

  std::map<uint64_t, BCSetting>           bc_settings;

public:
  explicit RadHydroSolver(const std::string& name);

  //01
  void Initialize() override;
  //02
  void Execute() override;
  //99
  void PrintRawOutput(const std::string& file_name);
  static VecDbl MakeOutputTimesFromStr(const std::string& output_times_str);

};

}//namespace chi_radhydro

#endif //RADHYDRO_SOLVER_BASE_H
