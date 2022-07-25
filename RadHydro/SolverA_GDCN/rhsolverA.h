#ifndef RADHYDRO_RHSOLVERA_H
#define RADHYDRO_RHSOLVERA_H

#include "../SolverBase/solver_base.h"

namespace chi_radhydro
{

class SolverA_GDCN : public RadHydroSolver
{
protected:
  double m_gamma;
  double m_Cv;

  std::string m_kappa_s_function;
  std::string m_kappa_a_function;

  std::vector<uint64_t> m_bndry_cells_local_ids;

  FVector m_material_advection;
  double m_radiation_advection;
public:
  //00
  explicit SolverA_GDCN(const std::string& name);
  //01
  void Initialize() override;
  //02
  void Execute() override;
  //02a
  void Predictor(SimRefs& sim_refs,
                 const std::vector<double>&      kappa_a_n,
                 const std::vector<double>&      kappa_t_n,
                 double dt,
                 const std::vector<UVector>&     U_n,
                 const std::vector<GradUTensor>& grad_U_n,
                 std::vector<UVector>&           U_nph,

                 const std::vector<double>&            rad_E_n,
                 const std::vector<chi_mesh::Vector3>& grad_rad_E_n,
                 std::vector<double>&                  rad_E_nph);

  //02b
  void CorrectorHydroAndMom(SimRefs& sim_refs,
                            const std::vector<double>&      kappa_a_n,
                            const std::vector<double>&      kappa_t_n,
                            const std::vector<double>&      kappa_a_nph,
                            const std::vector<double>&      kappa_t_nph,
                            double dt,
                            const std::vector<UVector>&     U_n,
                            const std::vector<UVector>&     U_nph,
                            const std::vector<GradUTensor>& grad_U_nph,
                            std::vector<UVector>&           U_np1,

                            const std::vector<double>&            rad_E_n,
                            const std::vector<double>&            rad_E_nph,
                            const std::vector<chi_mesh::Vector3>& grad_rad_E_nph,
                            std::vector<double>&                  rad_E_np1);

  static void AssembleGeneralEnergySystem(
    SimRefs&                        sim_refs,
    const std::vector<double>&      kappa_a_n,
    const std::vector<double>&      kappa_t_n,
    const std::vector<double>&      kappa_a_nph,
    const std::vector<double>&      kappa_t_nph,
    double                          tau,
    double                          theta1,
    double                          theta2,
    const std::vector<UVector>&     U_n,
    const std::vector<UVector>&     U_nph,
    const std::vector<UVector>&     U_nphstar,
    const std::vector<UVector>&     U_np1,
    const std::vector<GradUTensor>& grad_U_nph,
    const std::vector<double>&      rad_E_n,
    const std::vector<double>&      rad_E_nph,
    const std::vector<double>&      rad_E_nphstar,
    const std::vector<Vec3>&        grad_rad_E_nph,
    std::vector<double>&            k5_vec,
    std::vector<double>&            k6_vec,
    MatDbl &A, VecDbl &b);

  static void AssembleGeneralEnergySystem2(
    SimRefs&                        sim_refs,
    const std::vector<double>&      kappa_a_n,
    const std::vector<double>&      kappa_t_n,
    const std::vector<double>&      kappa_a_nph,
    const std::vector<double>&      kappa_t_nph,
    double                          tau,
    double                          theta1,
    double                          theta2,
    const std::vector<UVector>&     U_n,
    const std::vector<UVector>&     U_nph,
    const std::vector<UVector>&     U_nphstar,
    const std::vector<UVector>&     U_np1,
    const std::vector<GradUTensor>& grad_U_nph,
    const std::vector<double>&      rad_E_n,
    const std::vector<double>&      rad_E_nph,
    const std::vector<double>&      rad_E_nphstar,
    const std::vector<Vec3>&        grad_rad_E_nph,
    MatDbl &A, VecDbl &b);

  //02bb
  static double ComputeGradDotJ(SimRefs&                       sim_refs,
                                const chi_mesh::Cell&          cell_c,
                                      double                   sigma_t_c_n,
                                const std::vector<double>&     kappa_t_n,
                                const std::vector<UVector>&    U_n,
                                const std::vector<double>&     rad_E_n);

  static double MakeEmAbsSource(double sigma_a,
                         double T,
                         double radE);

  static double Make3rdGradRadEDot_u(
    SimRefs&                        sim_refs,
    const chi_mesh::Cell&           cell,
    double                          V_c,
    const std::vector<double>&      face_areas,
    const std::vector<double>&      rad_E,
    const std::vector<Vec3>&        grad_rad_E,
    const std::vector<UVector>&     U,
    const std::vector<GradUTensor>& grad_U);

  //91_utils
  void ProcessAndExportFields(const std::string& file_name,
                              const std::vector<UVector>& U,
                              const std::vector<double>& rad_E,
                              size_t num_local_nodes,
                              const std::map<std::string,double>& scalar_properties);

  struct SystemEnergy
  {
    double Emat = 0.0;
    double Erad = 0.0;
    double me_adv = 0.0;
    double re_adv = 0.0;
  };
  SystemEnergy ComputeSysEnergyChange(
                                double dt,
                                const std::vector<UVector>& U,
                                const std::vector<double>& rad_E,
                                const std::vector<GradUTensor>& grad_U,
                                const std::vector<chi_mesh::Vector3>& grad_rad_E,
                                const std::map<uint64_t, BCSetting>& bc_setttings);
};

}//namespace chi_radhydro

#endif //RADHYDRO_RHSOLVERA_H
