#ifndef RADHYDRO_SWEEPCHUNK_H
#define RADHYDRO_SWEEPCHUNK_H

namespace chi_radhydro
{
  class AngleSet;
}

namespace chi_radhydro
{
  class SweepChunk
  {
  public:
    virtual void Sweep(AngleSet& angle_set);
  };
}//namespace chi_radtran

#endif //RADHYDRO_SWEEPCHUNK_H
