#include "sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Initialize First-In-First-Out sweep scheduler.*/
void chi_radtran::SweepSchedulerFIFO::Initialize()
{
  chi::log.Log() << "Initializing SweepSchedulerFIFO";
}

//###################################################################
/**Execute First-In-First-Out sweep scheduler.*/
void chi_radtran::SweepSchedulerFIFO::Sweep(SweepChunk& chunk)
{
  for (auto& angle_set : m_angle_sets)
    chunk.Sweep(angle_set);
}