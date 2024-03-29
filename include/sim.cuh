#ifndef _SIM_H
#define _SIM_H

#include "parameter.h"
#include "particle.cuh"
#include "cell.cuh"
#include "data_pointer.cuh"

#include "thrust/device_vector.h"
#include "thrust/sort.h"
#include "glm/glm.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class Sim
{
public:
    Fluid_List fluid_list;
    Solid_List solid_list;
    Virt_List virt_list;
    thrust::device_vector<Cell> cell_list;

    Fluid_Host_List fluid_host_list;
    Solid_Host_List solid_host_list;
    Virt_Host_List virt_host_list;
    thrust::host_vector<Cell> cell_host_list;
    static Parameter parameter;
    static Data_Pointer data_pointer;

    void initiate()
    {
        // read file
        readParticleData();
        createCell();
        fillCell();
        copyData();
        setPointer();
        clearHost();
        std::cout << "\nInitialization completed" << std::endl;
    }

    void showMemoryUsed();

    void readParticleData();
    bool readFluidData(const std::string &path);
    bool readSolidData(const std::string &path);
    bool readVirtData(const std::string &path);

    void createCell();
    void fillCell();
    void copyData();
    int getCellID(const glm::vec3 &position);
    void setPointer();
    void clearHost();
};
#endif