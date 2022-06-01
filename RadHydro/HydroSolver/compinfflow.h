#ifndef RADHYDRO_COMPINFFLOW_H
#define RADHYDRO_COMPINFFLOW_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_matrixNXxNX.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_hydro
{

enum class CoordinateSystem
{
  NO_GEOMETRY_SET  = 0,
  ONED_SLAB        = 1,
  TWOD_CARTESIAN   = 4,
  THREED_CARTESIAN = 6
};

//Compressible Inviscid Fluid-Flow solver
//Comp         In       F     F ow
//CompInFFlow
class CompInFFlow : public chi_physics::Solver
{
protected:
  typedef std::shared_ptr<chi_math::SpatialDiscretization_FV> SDMFVPtr;
  typedef std::shared_ptr<chi_mesh::LogicalVolume> LogicalVolumePtr;
  typedef std::tuple<LogicalVolumePtr, std::string, double> LVFieldSetting;
  typedef chi_math::VectorN<5> UVector;
  typedef UVector FVector;
  typedef std::vector<UVector> GradUTensor;
  typedef chi_mesh::Vector3 Vec3;
  typedef chi_math::MatrixNXxNX<5,double> TMatrix;

public:
  double       t_max = -1.0;
  double       delta_t_max   = 1.0e-8;
  double       C_cfl         = 0.9;
  int64_t      num_timesteps = 1000;
  CoordinateSystem coordinate_system = CoordinateSystem::NO_GEOMETRY_SET;

  std::vector<LVFieldSetting> logvol_field_settings;

  chi_mesh::MeshContinuumPtr grid;
  SDMFVPtr                   fv;

  uint64_t num_nodes_local = 0;

  std::vector<UVector> U_old;

  std::vector<double> gamma, Cv;
  std::vector<double> rho, u, v, w, p, e;
  std::vector<double> temperature;
  std::vector<double> cell_char_length;

public:
  explicit CompInFFlow(const std::string& text_name);

  CompInFFlow(const CompInFFlow&) = delete;
  CompInFFlow operator=(const CompInFFlow&) = delete;

  void Initialize() override;
  //02
  void Execute() override;
  //02a
  double ComputeCourantLimitDelta_t() const;
  //02b
  std::vector<GradUTensor> ComputeGradients(const std::vector<UVector>& U) const;

  //90 utils
  static double PressureFromCellU(const UVector& U, double gamma);
  static FVector MakeF(const UVector& U,
                       double pressure,
                       const Vec3& n_f=Vec3(1,0,0));
  static double SoundSpeed(double gamma, double p, double rho);
  void FieldsToU(std::vector<UVector>& pU);
  void UToFields(const std::vector<UVector>& pU);

  static MatDbl MakeRotationMatrix(const Vec3& axis, double theta);
  static TMatrix MakeTransformationMatrix(const Vec3& n);
  static double MinMod(const std::vector<double>& a, bool verbose=false);
  static UVector MinModU(const std::vector<UVector>& vec_of_U,
                         bool verbose=false) ;
  static UVector UplusDXGradU(const UVector& U,
                              const Vec3& dx,
                              const GradUTensor& grad_U);

  static std::string PrintU(const UVector& U,
                            double epsilon=1.0e-12);

  //91
  static UVector HLLC_RiemannSolve(const UVector& U_L,
                                   const UVector& U_R,
                                   double p_L,
                                   double p_R,
                                   double gamma_L,
                                   double gamma_R,
                                   const Vec3& n_f,
                                   bool verbose=false) ;
};

}//namespace chi_radhydro

#endif //RADHYDRO_COMPINFFLOW_H
