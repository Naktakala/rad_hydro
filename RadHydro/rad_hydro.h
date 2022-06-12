#ifndef CHI_RADHYDRO_H
#define CHI_RADHYDRO_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_matrixNXxNX.h"

namespace chi_radhydro 
{
  static const double planck_constant_h = 6.62607015e-34;   // Joule/Hz
  static const double boltzmann_constant_kb = 1.380649e-23; // Joule/Kelvin
  static const double speed_of_light_mps =  299792458.0;
  static const double speed_of_light_cmpsh = 299.792;

  static const double a = 8.0*pow(M_PI,5)*pow(boltzmann_constant_kb,4)/
    (15.0*pow(planck_constant_h,3)*pow(speed_of_light_mps,3));

  typedef chi_math::VectorN<5> UVector;
  typedef chi_math::VectorN<5> FVector;
  typedef std::vector<UVector> GradUTensor;
  typedef chi_mesh::Vector3 Vec3;
  typedef chi_math::MatrixNXxNX<5,double> TMatrix;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;

  void TestFunction();

  void FieldsToU(std::vector<UVector> &pU,
                 const chi_mesh::MeshContinuum& grid,
                 const std::vector<double>& rho,
                 const std::vector<double>& u,
                 const std::vector<double>& v,
                 const std::vector<double>& w,
                 const std::vector<double>& e);

  double
  ComputeCourantLimitDelta_t(const chi_mesh::MeshContinuum& grid,
                             const std::vector<double>& rho,
                             const std::vector<double>& u,
                             const std::vector<double>& v,
                             const std::vector<double>& w,
                             const std::vector<double>& p,
                             const std::vector<double>& gamma,
                             const std::vector<double>& cell_char_length,
                             double C_cfl);
  double SoundSpeed(double gamma, double p, double rho);

  UVector MinModU(const std::vector<UVector> &vec_of_U,
                                bool verbose=false);

  double MinMod(const std::vector<double> &a,
                bool verbose=false);

  UVector
  UplusDXGradU(const UVector &U,
               const Vec3 &dx,
               const GradUTensor &grad_U);

  FVector
  MakeF(const UVector& U, double pressure, const Vec3& n_f=Vec3(1,0,0));

  chi_math::MatrixNXxNX<5,double> MakeTransformationMatrix(const Vec3 &n);

  MatDbl MakeRotationMatrix(const Vec3 &axis, double theta);
  double IdealGasPressureFromCellU(const UVector& U, double gamma);
  double IdealGasTemperatureFromCellU(const UVector& U, double C_v);
  double InternalEnergyFromCellU(const UVector& U);
  Vec3 VelocityFromCellU(const UVector& U);

  UVector
  HLLC_RiemannSolve(const UVector &U_L_raw,
                    const UVector &U_R_raw,
                    double p_L,
                    double p_R,
                    double gamma_L,
                    double gamma_R,
                    const Vec3& n_f,
                    bool verbose=false);

  std::string PrintU(const UVector &U, double epsilon=1.0e-12);

  void
  UToFields(const std::vector<UVector> &pU,
            const chi_mesh::MeshContinuum& grid,
            std::vector<double>& rho,
            std::vector<double>& u,
            std::vector<double>& v,
            std::vector<double>& w,
            std::vector<double>& e,
            std::vector<double>& p,
            const std::vector<double>& gamma);
}//namespace chi_radhydro 

#endif //CHI_RADHYDRO_H