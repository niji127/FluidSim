#ifndef _CELL_H
#define _CELL_H

#include "parameter.h"
#include "particle.cuh"
#include "global.h"

#include "thrust/device_vector.h"
#include "thrust/sort.h"

class Cell
{
public:
    int neighbor_cell_id[27];

    int fluid_start_id;
    int solid_start_id;
    int virt_start_id;

    int fluid_start_id_previous;
    int solid_start_id_previous;
    int virt_start_id_previous;

    size_t getMemoryUsed()
    {
        return sizeof(Cell);
    }

    Cell()
    {
        fluid_start_id = -1;
        solid_start_id = -1;
        virt_start_id = -1;

        fluid_start_id_previous = -1;
        solid_start_id_previous = -1;
        virt_start_id_previous = -1;

        for (int i = 0; i < 27;i++)
            neighbor_cell_id[i] = -1;
    }
};

#endif
