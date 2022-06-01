#ifndef RADHYDRO_SWEEPSCHEDULER_H
#define RADHYDRO_SWEEPSCHEDULER_H

#include "sweepchunk.h"
#include "angle_set.h"

#include <vector>

namespace chi_radtran
{
  class SweepScheduler
  {
  public:
    virtual void Initialize() = 0;
    virtual void Sweep(SweepChunk& chunk) = 0;
  };

  class SweepSchedulerFIFO : public SweepScheduler
  {
  private:
    std::vector<AngleSet>& m_angle_sets;
  public:
    explicit
    SweepSchedulerFIFO(std::vector<AngleSet>& angle_sets) :
    m_angle_sets(angle_sets) {}
    void Initialize() override;
    void Sweep(SweepChunk& chunk) override;
  };
}//namespace chi_radtran

#endif //RADHYDRO_SWEEPSCHEDULER_H
