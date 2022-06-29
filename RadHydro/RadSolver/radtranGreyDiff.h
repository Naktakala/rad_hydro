#ifndef RADHYDRO_RADTRANGREYDIFFUSION_H
#define RADHYDRO_RADTRANGREYDIFFUSION_H

#include "../rad_hydro.h"

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMath/chi_math.h"
#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_matrixNXxNX.h"

#include "ChiMesh/chi_mesh.h"

namespace chi_radhydro
{

class RadTranGreyDiffusion : public chi_physics::Solver
{
protected:
  typedef chi_math::SpatialDiscretization_FV SDM_FV;
  typedef std::shared_ptr<chi_mesh::LogicalVolume> LogicalVolumePtr;
  typedef std::tuple<LogicalVolumePtr, std::string, double> LVFieldSetting;

  std::shared_ptr<chi_mesh::MeshContinuum> grid;

private:
  std::shared_ptr<SDM_FV>                  fv;
  chi_math::UnknownManager                 unitary_uk_man;

public:
  std::vector<LVFieldSetting> logvol_field_settings;

  double       t_max         = -1.0;
  double       delta_t_max   = 1.0e-8;
  double       C_cfl         = 0.9;
  int64_t      num_timesteps = 1000;


  uint64_t num_nodes_local = 0;

  std::vector<UVector> U_old;

  std::vector<double> gamma, Cv;
  std::vector<double> rho, u, v, w, p, e, rad_E;
  std::vector<double> temperature;
  std::vector<double> cell_char_length;

public:
  //00 Constructor
  explicit RadTranGreyDiffusion(const std::string& name);

  RadTranGreyDiffusion(const RadTranGreyDiffusion&) = delete;
  RadTranGreyDiffusion operator=(const RadTranGreyDiffusion&) = delete;

  //01
  void Initialize() override;
  //02
  void Execute() override;
  void Predictor(double dt,
                 std::vector<UVector>&           U_n,
                 std::vector<UVector>&           U_nph,
                 std::vector<GradUTensor>&       grad_U,

                 std::vector<double>&            rad_E_n,
                 std::vector<double>&            rad_E_nph,
                 std::vector<chi_mesh::Vector3>& grad_rad_E);
  void Corrector(double dt,
                 std::vector<UVector>&           U_n,
                 std::vector<UVector>&           U_nph,
                 std::vector<UVector>&           U_np1,
                 std::vector<GradUTensor>&       grad_U,

                 std::vector<double>&            rad_E_n,
                 std::vector<double>&            rad_E_nph,
                 std::vector<double>&            rad_E_np1,
                 std::vector<chi_mesh::Vector3>& grad_rad_E);
  //02a
  static void AssemblePredictorRadESystem(MatDbl &A, VecDbl &b,
                                          chi_mesh::MeshContinuum&  grid_ref,
                                          std::shared_ptr<SDM_FV>&  fv_ref,
                                          std::vector<UVector>&     U_n,
                                          std::vector<UVector>&     U_nstar,
                                          std::vector<UVector>&     U_nph,
                                          std::vector<double>&      rad_E_n,
                                          std::vector<double>&      rad_E_nstar,
                                          std::vector<double>&      kappa_t,
                                          std::vector<double>&      kappa_a,
                                          std::vector<GradUTensor>& grad_U_n,
                                          std::vector<Vec3>&        grad_rad_E,
                                          std::vector<double>&      Cv,
                                          std::vector<double>&      k5_vec,
                                          std::vector<double>&      k6_vec,
                                          double                    delta_t);

  //03
  std::vector<GradUTensor>
  ComputeUGradients(const std::vector<UVector>& U) const;

  std::vector<Vec3> ComputeRadEGradients(const std::vector<double>& E) const;
};

}//namespace chi_radtran

#endif //RADHYDRO_RADTRANGD_H
