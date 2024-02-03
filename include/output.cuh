#ifndef _OUTPUT_CUH
#define _OUTPUT_CUH

#include "function.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Output : public DeviceFunction
{
public:
    static void writeFluidVelocity(Sim *sim, const std::string &path, SimTime *time);
    static void writeFluidDensity(Sim *sim, const std::string &path, SimTime *time);
    static void writeFluidPressure(Sim *sim, const std::string &path, SimTime *time);

    static void writeSolidVelocity(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidPressure(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidDensity(Sim *sim, const std::string &path, SimTime *time);
    static void writeSolidNormal(Sim *sim, const std::string &path, SimTime *time);

    static void writeVirtNormal(Sim *sim, const std::string &path, SimTime *time);
    static void writeVirtPressure(Sim *sim, const std::string &path, SimTime *time);

    static void outputData(Sim *sim);
};

#endif