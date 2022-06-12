#ifndef RADTRANDO_H
#define RADTRANDO_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMath/chi_math.h"
#include "ChiMesh/chi_mesh.h"

#include "SPDS/sweepstructure.h"
#include "SPDS/angle_set.h"
#include "SPDS/sweepscheduler.h"

namespace chi_radhydro
{

class RadTranDO : public chi_physics::Solver
{
private:
  typedef chi_math::SpatialDiscretization_PWLD SDM_PWLD;
private:
  struct
  {
    size_t num_angles = 8;
    size_t scattering_order = 0;
    size_t num_groups = 0;
  }options;

  std::shared_ptr<chi_math::AngularQuadrature> angular_quadrature;
  std::shared_ptr<chi_mesh::MeshContinuum> grid;
  chi_math::UnknownManager uk_man_phi;
  std::shared_ptr<SDM_PWLD> sdm;

  size_t num_moments = 1;
  uint64_t num_local_phiDOFs = 0;

  std::vector<double> phi_old;

  std::vector<AngleSet> angle_sets;
  std::unique_ptr<SweepScheduler> sweep_scheduler;

public:
  //00 Constructor
  explicit RadTranDO(const std::string& name);

  RadTranDO(const RadTranDO&) = delete;
  RadTranDO operator=(const RadTranDO&) = delete;

  //01 Initialize
  void Initialize() override;

  //02 Execute
  void Execute() override;
};
}//namespace chi_radtran


#endif// RADTRANDO_H

