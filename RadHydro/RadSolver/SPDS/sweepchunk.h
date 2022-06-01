#ifndef RADHYDRO_SWEEPCHUNK_H
#define RADHYDRO_SWEEPCHUNK_H

namespace chi_radtran
{
  class AngleSet;
}

namespace chi_radtran
{
  class SweepChunk
  {
  public:
    virtual void Sweep(AngleSet& angle_set);
  };
}//namespace chi_radtran

#endif //RADHYDRO_SWEEPCHUNK_H
