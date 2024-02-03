#include "function.cuh"
#include "calculation.cuh"
#include "memory_copy.cuh"
#include "sim.cuh"
#include "output.cuh"
#include "probe.h"

#include "thrust/host_vector.h"

#include <fstream>
#include <string>
#include <windows.h>

void DeviceFunction::GPUCalulation(Sim *sim)
{
    MemoryCopy::copyConstant();
    DeviceCalculation::calculation(sim);
}

void DeviceCalculation::calculation(Sim *sim)
{
    SimTime *time = Sim::parameter.getTime();
    while (time->addStep())
    {
        Output::outputData(sim);
        particleCellUpdate();
        // sortParticle(sim);
        fillCell();
        coupleCalculation();
        virtUpdate();
        fluidCalculation();
        solidCalculation(sim);
    }
    Output::outputData(sim);
}