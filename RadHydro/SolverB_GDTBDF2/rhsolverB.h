#ifndef RADHYDRO_RHSOLVERB_H
#define RADHYDRO_RHSOLVERB_H

#include "../SolverA_GDCN/rhsolverA.h"

namespace chi_radhydro
{

class SolverB_GDTBDF : public SolverA_GDCN
{
protected:
public:
  //00
  explicit SolverB_GDTBDF(const std::string& name);
  //02
  void Execute() override;
  //02b
  void CorrectorB(const std::map<uint64_t, BCSetting>& bc_setttings,
                  const std::vector<double>&      kappa_a_n,
                  const std::vector<double>&      kappa_t_n,
                  const std::vector<double>&      kappa_a_nph,
                  const std::vector<double>&      kappa_t_nph,
                  const std::vector<double>&      kappa_a_np3q,
                  const std::vector<double>&      kappa_t_np3q,
                  double dt,
                  const std::vector<UVector>&     U_n,
                  const std::vector<UVector>&     U_nph,
                  const std::vector<UVector>&     U_np3q,
                  const std::vector<GradUTensor>& grad_U_np3q,
                  std::vector<UVector>&           U_np1,

                  const std::vector<double>&            rad_E_n,
                  const std::vector<double>&            rad_E_nph,
                  const std::vector<double>&            rad_E_np3q,
                  const std::vector<chi_mesh::Vector3>& grad_rad_E_np3q,
                  std::vector<double>&                  rad_E_np1);
  //02ba
  static void AssembleGeneralEnergySystemB(
    const chi_mesh::MeshContinuum&  grid_ref,
    std::shared_ptr<SDM_FV>&        fv_ref,
    const std::map<uint64_t, BCSetting>& bc_setttings,
    const std::vector<double>&      kappa_a_n,
    const std::vector<double>&      kappa_t_n,
    const std::vector<double>&      kappa_a_nph,
    const std::vector<double>&      kappa_t_nph,
    const std::vector<double>&      kappa_a_np3q,
    const std::vector<double>&      kappa_t_np3q,
    double                          Cv,
    double                          tau,
    double                          theta1,
    double                          theta2,
    double                          theta3,

    const std::vector<UVector>&     U_n,
    const std::vector<UVector>&     U_nph,
    const std::vector<UVector>&     U_np3q,
    const std::vector<UVector>&     U_np3qstar,
    const std::vector<UVector>&     U_np1,
    const std::vector<GradUTensor>& grad_U_np3q,

    const std::vector<double>&      rad_E_n,
    const std::vector<double>&      rad_E_nph,
    const std::vector<double>&      rad_E_np3q,
    const std::vector<double>&      rad_E_np3qstar,
    const std::vector<Vec3>&        grad_rad_E_np3q,
    std::vector<double>&            k5_vec,
    std::vector<double>&            k6_vec,
    MatDbl &A, VecDbl &b);
};

}//namespace chi_radhydro

#endif //RADHYDRO_RHSOLVERB_H
