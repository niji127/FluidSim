#include "calculation.cuh"
#include "memory_copy.cuh"
#include "global.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ int findPositionCell(const glm::vec3 &pos)
{
    float xmin = dev_par.domain.cell_min[0];
    float ymin = dev_par.domain.cell_min[1];
    float zmin = dev_par.domain.cell_min[2];
    float inter_inv_x = dev_par.domain.interval_inv[0];
    float inter_inv_y = dev_par.domain.interval_inv[1];
    float inter_inv_z = dev_par.domain.interval_inv[2];
    int block_ny = dev_par.domain.cell_number[1];
    int block_nz = dev_par.domain.cell_number[2];

    int nx = (int)((pos.x - xmin) * inter_inv_x);
    int ny = (int)((pos.y - ymin) * inter_inv_y);
    int nz = (int)((pos.z - zmin) * inter_inv_z);

    return nx * block_ny * block_nz + ny * block_nz + nz;
}

__global__ void fillFluidCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int fluid_num = dev_par.fluid.number;
    if (id >= fluid_num)
        return;
    dev_data.fluid.cell_id[id] = findPositionCell(dev_data.fluid.position[id]);
}

__global__ void fillSolidCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int solid_num = dev_par.solid.number;
    if (id >= solid_num)
        return;
    dev_data.solid.cell_id[id] = findPositionCell(dev_data.solid.position[id]);
}

__global__ void fillVirtCellID()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int virt_num = dev_par.virt.number;
    if (id >= virt_num)
        return;
    dev_data.virt.cell_id[id] = findPositionCell(dev_data.virt.position[id]);
}

void DeviceCalculation::particleCellUpdate()
{
    int thread_num = Sim::parameter.kernel.thread_num;

    int fluid_num = Sim::parameter.fluid.number;
    int fluid_block = fluid_num / thread_num + 1;
    fillFluidCellID<<<fluid_block, thread_num>>>();

    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;
    fillSolidCellID<<<solid_block, thread_num>>>();

    int virt_num = Sim::parameter.virt.number;
    int virt_block = virt_num / thread_num + 1;
    fillVirtCellID<<<virt_block, thread_num>>>();

    cudaDeviceSynchronize();
}

__global__ void fillCellList()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.domain.cell_number_total)
        return;

    Cell *cell_i = dev_data.cell + id;
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virt_Data virt = dev_data.virt;

    int current_fluid = -1;
    int current_solid = -1;
    int current_virt = -1;
    
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int fluid_id = cell_j->fluid_start_id_previous;
        while (fluid_id != -1)
        {
            if (fluid.cell_id[fluid_id] == id)
            {
                if (cell_i->fluid_start_id == -1)
                {
                    cell_i->fluid_start_id = fluid_id;
                    current_fluid = fluid_id;
                }
                else
                {
                    fluid.next_id[current_fluid] = fluid_id;
                    current_fluid = fluid_id;
                }
            }
            fluid_id = fluid.next_id_previous[fluid_id];
        }

        int solid_id = cell_j->solid_start_id_previous;
        while (solid_id != -1)
        {
            if (solid.cell_id[solid_id] == id)
            {
                if (cell_i->solid_start_id == -1)
                {
                    cell_i->solid_start_id = solid_id;
                    current_solid = solid_id;
                }
                else
                {
                    solid.next_id[current_solid] = solid_id;
                    current_solid = solid_id;
                }
            }
            solid_id = solid.next_id_previous[solid_id];
        }

        int virt_id = cell_j->virt_start_id_previous;
        while (virt_id != -1)
        {
            if (virt.cell_id[virt_id] == id)
            {
                if (cell_i->virt_start_id == -1)
                {
                    cell_i->virt_start_id = virt_id;
                    current_virt = virt_id;
                }
                else
                {
                    virt.next_id[current_virt] = virt_id;
                    current_virt = virt_id;
                }
            }
            virt_id = virt.next_id_previous[virt_id];
        }
    }
}
__global__ void copyCellList()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.domain.cell_number_total)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virt_Data virt = dev_data.virt;

    Cell *cell = dev_data.cell + id;

    int fluid_id = cell->fluid_start_id;
    while (fluid_id != -1)
    {
        fluid.next_id_previous[fluid_id] = fluid.next_id[fluid_id];
        fluid.next_id[fluid_id] = -1;
        fluid_id = fluid.next_id_previous[fluid_id];
    }
    cell->fluid_start_id_previous = cell->fluid_start_id;
    cell->fluid_start_id = -1;

    int solid_id = cell->solid_start_id;
    while (solid_id != -1)
    {
        solid.next_id_previous[solid_id] = solid.next_id[solid_id];
        solid.next_id[solid_id] = -1;
        solid_id = solid.next_id_previous[solid_id];
    }
    cell->solid_start_id_previous = cell->solid_start_id;
    cell->solid_start_id = -1;

    int virt_id = cell->virt_start_id;
    while (virt_id != -1)
    {
        virt.next_id_previous[virt_id] = virt.next_id[virt_id];
        virt.next_id[virt_id] = -1;
        virt_id = virt.next_id_previous[virt_id];
    }
    cell->virt_start_id_previous = cell->virt_start_id;
    cell->virt_start_id = -1;
}

void DeviceCalculation::fillCell()
{
    int thread_num = Sim::parameter.kernel.thread_num;
    int cell_num = Sim::parameter.domain.cell_number_total;
    int cell_block = cell_num / thread_num + 1;
    SimTime *time = Sim::parameter.getTime();

    copyCellList<<<cell_block, thread_num>>>();
    cudaDeviceSynchronize();

    fillCellList<<<cell_block, thread_num>>>();
    cudaDeviceSynchronize();
}