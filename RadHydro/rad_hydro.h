#ifndef CHI_RADHYDRO_H
#define CHI_RADHYDRO_H

#include "ChiMesh/chi_mesh.h"

#include "ChiMath/chi_math.h"
#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_matrixNXxNX.h"

#include "ChiTimer/chi_timer.h"

#include <map>



namespace chi_radhydro 
{
  static const double boltzmann_constant_kb  = 1.380649e-23;   // Joule/Kelvin
  static const double boltzmann_constant_kb1 = 1.60223860E-25; // Jk/keV
  static const double boltzmann_constant_kb2 = 8.617e-8; // keV/Kelvin

  static const double planck_constant_h  = 6.62607015e-34;   // Joule/Hz
  static const double planck_constant_h1 = 6.62607015e-35;   // Jk-sh
  static const double planck_constant_h2 = 4.13380133E-10;   // keV-sh

//  static const double speed_of_light_mps =  299792458.0;
  static const double speed_of_light_cmpsh = 299.792458;

//  static const double a = (8.0/15.0)*pow(M_PI,5)*pow(boltzmann_constant_kb1,4)/
//                          (pow(planck_constant_h1*speed_of_light_cmpsh,3));
//  static const double a = 0.013722354852;
  static const double a = 0.0137201720;

  static const chi_mesh::Vector3 ihat(1,0,0);
  static const chi_mesh::Vector3 jhat(0,1,0);
  static const chi_mesh::Vector3 khat(0,0,1);

  enum UVectorEntries
  {
    RHO   = 0,
    RHO_U = 1,
    RHO_V = 2,
    RHO_W = 3,
    MAT_E = 4
  };

  enum class BCType : int
  {
    TRANSMISSIVE = 0,
    FIXED        = 1
  };

  struct BCSetting
  {
    BCType                type = BCType::TRANSMISSIVE;
    std::array<double, 7> values = {0,0,0,0,0,0,0};
  };

  typedef chi_math::VectorN<5> UVector;
  typedef chi_math::VectorN<5> FVector;
  typedef std::vector<UVector> GradUTensor;
  typedef chi_mesh::Vector3 Vec3;
  typedef chi_math::MatrixNXxNX<5,double> TMatrix;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  typedef std::pair<TMatrix, TMatrix> TMatrixPair;

  struct SimRefs
  {
    const chi_mesh::MeshContinuum&       grid;
    chi_math::SpatialDiscretization_FV&  fv;
    const std::map<uint64_t, BCSetting>& bc_settings;
    const double                         Cv;
    const double                         gamma;
    const std::string&                   kappa_s_function;
    const std::string&                   kappa_a_function;
  };

  //00_general
  double SoundSpeed(double gamma, double p, double rho);

  double ComputeCourantLimitDelta_t(const std::vector<UVector>& vecU,
                                    double gamma,
                                    const std::vector<double>& cell_char_length,
                                    double C_cfl);

  UVector MakeUFromBC(const BCSetting& bc_setting,
                      const UVector& U_default);
  double MakeRadEFromBC(const BCSetting& bc_setting,
                        double radE_default);
  double GetResetTimer(chi_objects::ChiTimer& timer);


  //01_fieldsToFrom
  void FieldsToU(const std::vector<double>& rho,
                 const std::vector<double>& u,
                 const std::vector<double>& v,
                 const std::vector<double>& w,
                 const std::vector<double>& e,
                 std::vector<UVector> &pU);

  void
  UToFields(const std::vector<UVector> &pU,
            std::vector<double>& rho,
            std::vector<double>& u,
            std::vector<double>& v,
            std::vector<double>& w,
            std::vector<double>& e,
            std::vector<double>& p,
            double gamma);

  //02_Ustuff
  UVector
    UplusDXGradU(const UVector &U, const Vec3 &dx,const GradUTensor &grad_U);
  double IdealGasPressureFromCellU(const UVector& U, double gamma);
  double IdealGasTemperatureFromCellU(const UVector& U, double C_v);
  double InternalEnergyFromCellU(const UVector& U);
  Vec3 VelocityFromCellU(const UVector& U);
  std::string PrintU(const UVector &U, double epsilon=1.0e-12);

  //03_Fstuff
  MatDbl MakeRotationMatrix(const Vec3 &axis, double theta);
  chi_math::MatrixNXxNX<5,double> MakeTransformationMatrix(const Vec3 &n);
  TMatrixPair MakeTandTinv(const Vec3& n);
  FVector MakeF(const UVector& U, double pressure, const Vec3& n_f=Vec3(0,0,1));
  FVector MakeFNoTransform(const UVector& U, double pressure);
  FVector Fswap_w_and_u(const FVector &F);
  UVector Uswap_w_and_u(const UVector &U);
  FVector MakeFWithRadE(const UVector& U, double pressure,
                        double radE, const Vec3& n_f=Vec3(0,0,1));


  //04_gradients
  UVector MinModU(const std::vector<UVector> &vec_of_U,
                                bool verbose=false);

  double MinMod(const std::vector<double> &ais,
                bool verbose=false);

  std::vector<GradUTensor>
    ComputeUGradients(const std::vector<UVector>& U, SimRefs& sim_refs);

  std::vector<Vec3>
    ComputeRadEGradients(const std::vector<double>& radE, SimRefs& sim_refs);

  //05_Riemann
  UVector
  HLLC_RiemannSolve(const UVector &U_L_raw,
                    const UVector &U_R_raw,
                    double gamma,
                    const Vec3& n_f,
                    bool verbose=false);

  //06_kappa
  double ComputeKappaFromLua(double T, int mat_id, const std::string& lua_fname);

  void ComputeCellKappas(const chi_mesh::MeshContinuum& grid,
                         double                         Cv,
                         const std::vector<UVector>&    U,
                         const std::string&             kappa_s_function,
                         const std::string&             kappa_a_function,
                         std::vector<double>&     kappa_a,
                         std::vector<double>&     kappa_t);

  //07_MHM
  void MHM_HydroRadEPredictor(
    SimRefs&                              sim_refs,
    double                                tau,
    const std::vector<UVector>&           U_a,
    const std::vector<GradUTensor>&       grad_U_a,
    const std::vector<double>&            rad_E_a,
    const std::vector<chi_mesh::Vector3>& grad_rad_E_a,
    std::vector<UVector>&                 U_a_star,
    std::vector<double>&                  rad_E_a_star
    );

  void MHM_HydroRadECorrector(
    SimRefs&                              sim_refs,
    double                                tau,
    const std::vector<UVector>&           U_old,
    const std::vector<UVector>&           U_int,
    const std::vector<GradUTensor>&       grad_U_int,
    const std::vector<double>&            rad_E_old,
    const std::vector<double>&            rad_E_int,
    const std::vector<chi_mesh::Vector3>& grad_rad_E_int,
    std::vector<UVector>&                 U_int_star,
    std::vector<double>&                  rad_E_int_star
    );

  void DensityMomentumUpdateWithRadMom(
    SimRefs&                              sim_refs,
    const std::vector<double>&            kappa_t,
    double                                tau,
    const std::vector<UVector>&           U_old,
    const std::vector<UVector>&           U_int,
    const std::vector<double>&            rad_E_old,
    const std::string&                    kappa_s_function,
    const std::string&                    kappa_a_function,
    std::vector<UVector>&                 U_new
  );

  //08
  std::vector<double> MakePostShockConditionsRH(
    double Cv,
    double gamma,
    double rho0,
    double T0,
    double u0,
    double rho1_guess,
    double T1_guess,
    double u1_guess);
  std::vector<double> MakePostShockConditionsHydroOnly(
    double Cv,
    double gamma,
    double rho0,
    double T0,
    double u0,
    double rho1_guess,
    double T1_guess,
    double u1_guess);
}//namespace chi_radhydro 

#endif //CHI_RADHYDRO_H