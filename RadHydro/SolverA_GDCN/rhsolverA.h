#ifndef RADHYDRO_RHSOLVERA_H
#define RADHYDRO_RHSOLVERA_H

#include "../SolverBase/solver_base.h"

namespace chi_radhydro
{

class SolverA_GDCN : public RadHydroSolver
{
protected:
public:
  //00
  explicit SolverA_GDCN(const std::string& name);
  //01
  void Initialize() override;
  //02
  void Execute() override;
  //02a
  void Predictor(const std::map<int, BCSetting>& bc_setttings,
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
  void Corrector(const std::map<int, BCSetting>& bc_setttings,
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
    const chi_mesh::MeshContinuum&  grid_ref,
    std::shared_ptr<SDM_FV>&        fv_ref,
    const std::map<int, BCSetting>& bc_setttings,
    const std::vector<double>&      kappa_a_n,
    const std::vector<double>&      kappa_t_n,
    const std::vector<double>&      kappa_a_nph,
    const std::vector<double>&      kappa_t_nph,
    const std::vector<double>&      Cv,
    double                          tau,
    double                          theta1,
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
};

}//namespace chi_radhydro

#endif //RADHYDRO_RHSOLVERA_H
