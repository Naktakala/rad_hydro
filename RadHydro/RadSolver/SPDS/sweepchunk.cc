#include "sweepchunk.h"

#include "angle_set.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Execution of a sweep chunk.*/
void chi_radtran::SweepChunk::Sweep(AngleSet &angle_set)
{
  const auto& sweep_structure = angle_set.GetSweepStructure();
  const auto& sweep_ordering  = sweep_structure.GridOrdering();

  for (size_t cell_local_id : sweep_ordering)
  {
    chi::log.Log() << "Cell " << cell_local_id;
    for (auto& angle : angle_set.Angles())
    {
      chi::log.Log() << "Sweeping angle " << angle.id;
    }//for angle in angle-set
  }//for cell in sweep order
}