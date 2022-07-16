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
  std::shared_ptr<chi_mesh::MeshContinuum> m_grid;
  std::shared_ptr<SDM_FV>                  m_fv;

public:
  std::vector<LVFieldSetting>              m_logvol_field_settings;
  std::map<std::string, ScalarField>       m_scalar_fields;

  std::map<uint64_t, BCSetting>            m_bc_settings;

public:
  explicit RadHydroSolver(const std::string& name);

  //01
  void Initialize() override;
  //02
  void Execute() override;
  //99
  void PrintRawOutput(const std::string& file_name,
                      const std::map<std::string,double>& scalar_properties);
  static VecDbl MakeOutputTimesFromStr(const std::string& output_times_str);

};

}//namespace chi_radhydro

#endif //RADHYDRO_SOLVER_BASE_H
