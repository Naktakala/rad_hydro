#ifndef RADHYDRO_ANGLE_SET_H
#define RADHYDRO_ANGLE_SET_H

#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "sweepstructure.h"

namespace chi_radtran
{
  class AngularQuadraturePointInfo
  {
  public:
    chi_math::QuadraturePointPhiTheta phi_theta;
    chi_mesh::Vector3                 omega;
    double                            weight;
    size_t                            id;

    AngularQuadraturePointInfo(
      chi_math::QuadraturePointPhiTheta p_phi_theta,
      const chi_mesh::Vector3&          p_omega,
      double                            p_weight,
      size_t                            p_id) :
      phi_theta(p_phi_theta),
      omega(p_omega),
      weight(p_weight),
      id(p_id) {}
  };

  class AngleSet
  {
  private:
    std::vector<AngularQuadraturePointInfo> m_angles;
    std::shared_ptr<SweepStructure>         m_sweep_structure;
  public:
    explicit
    AngleSet(std::vector<AngularQuadraturePointInfo> angle_points_info,
             std::shared_ptr<SweepStructure> sweep_ordering) :
      m_angles(std::move(angle_points_info)),
      m_sweep_structure(std::move(sweep_ordering)){}

    const std::vector<AngularQuadraturePointInfo>& Angles() const
    {return m_angles;}
    const SweepStructure& GetSweepStructure() const {return *m_sweep_structure;}
  };
}//namespace chi_radtran

#endif //RADHYDRO_ANGLE_SET_H
