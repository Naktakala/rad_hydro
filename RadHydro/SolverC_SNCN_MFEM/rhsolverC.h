#ifndef RADHYDRO_RHSOLVERC_H
#define RADHYDRO_RHSOLVERC_H

#include "../SolverA_GDCN/rhsolverA.h"
#include "ChiMath/chi_math.h"

namespace chi_radhydro
{

class SolverC_SNCN_MFEM : public SolverA_GDCN
{
protected:
  typedef chi_math::SpatialDiscretization_PWLC SDM_PWLC;
  typedef chi_math::SpatialDiscretization_PWLD SDM_PWLD;
protected:
  std::shared_ptr<SDM_PWLC> m_pwlc;
  std::shared_ptr<SDM_PWLD> m_pwld;

  std::vector<MatDbl> m_list_of_Cc_star;
  std::vector<MatDbl> m_list_of_Cvol_star;
  size_t              m_num_local_cfem_nodes=0;
  size_t              m_num_local_dfem_nodes=0;
public:
  //00
  explicit SolverC_SNCN_MFEM(const std::string& name);
  //01
  void Initialize() override;
  //02
  void Execute() override;
    //02a
  void FieldsToRadEF(const std::vector<UVector>& U_n,
                     VecDbl& rad_E_n,
                     std::vector<Vec3>& rad_F_n,
                     std::vector<Vec3>& rad_F0_n);
    //02b
  void DensityMomentumUpdateWithRadF0(
      SimRefs&                              sim_refs,
      const std::vector<double>&            kappa_t,
      double                                tau,
      const std::vector<UVector>&           U_old,
      const std::vector<UVector>&           U_int,
      const std::vector<Vec3>&              rad_F0_old,
      std::vector<UVector>&                 U_new
    );
  //02ca
  static void AssySolveMFEMEnergyExchange(
    SimRefs&                        sim_refs,
    SDM_PWLC&                       pwl,
    SDM_PWLD&                       pwld,
    std::vector<MatDbl>&            list_of_Cc_star,
    std::vector<MatDbl>&            list_of_Cvol_star,
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
    const std::vector<Vec3>&        grad_rad_E_nph,
    const std::vector<Vec3>&        rad_F0_nph,
    const std::vector<Vec3>&        rad_F_n,
    const std::vector<double>&      VEFf_nodal,
    const std::vector<double>&      VEFf_ctr,
    std::vector<UVector>&           U_np1,
    std::vector<double>&            rad_E_np1,
    std::vector<Vec3>&              rad_F_np1,
    std::vector<Vec3>&              rad_F0_np1);

  //03
  static void Sweep1D(SimRefs&                           sim_refs,
                      SDM_PWLC&                          pwlc,
                      SDM_PWLD&                          pwld,
                      const std::vector<double>&         kappa_a,
                      const std::vector<double>&         kappa_t,
                      const std::vector<UVector>&        U,
                      const std::vector<double>&         radE,
                      const std::vector<Vec3>&           radF0,
                      const chi_math::AngularQuadrature& quadrature,
                      VecDbl& VEFf_nodal,
                      VecDbl& VEFf_ctr,
                      bool verbose = false);


  };

}//namespace chi_radhydro

#endif //RADHYDRO_RHSOLVERC_H
