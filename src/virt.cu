#include "calculation.cuh"
#include "memory_copy.cuh"
#include "global.h"
#include <glm/gtx/norm.hpp>

extern __device__ float densityToPressure(float &density);
extern __device__ float pressureToDensity(float &pressure);

__global__ void getVirtData()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.virt.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Virt_Data virt = dev_data.virt;
    Cell *cell_i = dev_data.cell + virt.cell_id[id];
    glm::vec3 gravity(dev_par.physics.gravity[0], dev_par.physics.gravity[1], dev_par.physics.gravity[2]);

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    float coef_diff = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 velocity = glm::vec3(0.0f);
    glm::vec3 normal = glm::vec3(0.0f);
    float pressure(0.0f), density(dev_par.fluid.reference_density);
    float nearest_distance_square(LARGE_FLOAT);

    glm::vec3 x_i = virt.position[id];
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int id_j = cell_j->fluid_start_id;
        while (id_j!=-1)
        {
            glm::vec3 dx = fluid.position[id_j] - x_i;
            float dis_2 = glm::length2(dx);
            if (dis_2 < nearest_distance_square)
            {
                density = fluid.density[id_j];
                pressure = fluid.pressure[id_j] - density * glm::dot(gravity, dx);
                velocity = fluid.velocity[id_j];
                nearest_distance_square = dis_2;
            }
            id_j = fluid.next_id[id_j];
        }

        id_j = cell_j->virt_start_id;
        while (id_j != -1)
        {
            glm::vec3 dx = virt.position[id_j] - x_i;
            float dis_2 = glm::length2(dx);
            if (dis_2 < impact_dis_2 && id != id_j)
            {
                float q_2 = dis_2 * h_inv_2;
                float kernel_coef_diff = expf(-q_2) * coef_diff;
                normal += -dx * kernel_coef_diff;
            }
            id_j = virt.next_id[id_j];
        }
    }

    if (glm::length2(normal) < EPS_FOR_SUM)
        normal = glm::vec3(0.0f);
    else
        normal = glm::normalize(normal);
    virt.normal[id] = normal;
    virt.pressure[id] = pressure;
    virt.velocity[id] = velocity - 2.0f * glm::dot(velocity, normal) * normal;
    virt.density[id] = density;
}

void DeviceCalculation::virtUpdate()
{
    if (!(Sim::parameter.hasVirtual()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int virt_num = Sim::parameter.virt.number;
    int virt_block = virt_num / thread_num + 1;
    getVirtData<<<virt_block, thread_num>>>();
    cudaDeviceSynchronize();
}