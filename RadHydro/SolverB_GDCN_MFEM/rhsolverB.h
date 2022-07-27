#ifndef RADHYDRO_RHSOLVERB_H
#define RADHYDRO_RHSOLVERB_H

#include "../SolverA_GDCN/rhsolverA.h"
#include "ChiMath/chi_math.h"

namespace chi_radhydro
{

class SolverB_GDCN_MFEM : public SolverA_GDCN
{
protected:
  typedef chi_math::SpatialDiscretization_PWLC SDM_PWLC;
protected:
  std::shared_ptr<SDM_PWLC> m_pwl;

  std::vector<MatDbl> m_list_of_Cc_star;
  size_t              m_num_local_cfem_nodes;
public:
  //00
  explicit SolverB_GDCN_MFEM(const std::string& name);
  //01
  void Initialize() override;
  //02
  void Execute() override;
  //02ba
  static void AssembleMFEMEnergySystem(
    SimRefs&                        sim_refs,
    SDM_PWLC&                       pwl,
    std::vector<MatDbl>&            list_of_Cc_star,
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

  static void AssySolveMFEMEnergyExchange(
    SimRefs&                        sim_refs,
    SDM_PWLC&                       pwl,
    std::vector<MatDbl>&            list_of_Cc_star,
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
    const std::vector<GradUTensor>& grad_U_nph,
    const std::vector<double>&      rad_E_n,
    const std::vector<double>&      rad_E_nph,
    const std::vector<double>&      rad_E_nphstar,
    const std::vector<Vec3>&        grad_rad_E_nph,
    std::vector<UVector>&           U_np1,
    std::vector<double>&            rad_E_np1);

};

}//namespace chi_radhydro

#endif //RADHYDRO_RHSOLVERB_H
