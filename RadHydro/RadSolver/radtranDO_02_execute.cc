#include "radtranDO.h"

void chi_radhydro::RadTranDO::Execute()
{
  //Lets do classic richardson first

  SweepChunk test_chunk;

  for (size_t k=0; k<1; ++k)
  {
    //Set RHS
    sweep_scheduler->Sweep(test_chunk);
  }//For iteration k
}