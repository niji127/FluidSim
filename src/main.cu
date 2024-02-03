// SPH solver 2023.6

#include "sim.cuh"
#include "function.cuh"
#include "cuda_macro.cuh"
#include "global.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <windows.h>

#include "thrust/iterator/constant_iterator.h"
#include "thrust/host_vector.h"

int main()
{
	// initialize
	Sim::parameter.initiate();
	Sim::parameter.showInfo();

	Sim *sim = new Sim;
	sim->initiate();
	sim->showMemoryUsed();

	DeviceFunction::GPUCalulation(sim);

	return 0;
}