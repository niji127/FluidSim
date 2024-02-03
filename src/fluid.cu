#include "calculation.cuh"
#include "memory_copy.cuh"

#include "thrust/host_vector.h"
#include <glm/gtx/norm.hpp>

__device__ float pressureToDensity(float &pressure)
{
    float density_reference = dev_par.fluid.reference_density;
    float coef = dev_par.fluid.coefficient_p2rho;
    float gamma_inv = dev_par.fluid.gamma_inv;

    return density_reference * powf(pressure * coef + 1.0f, gamma_inv);
}

__device__ float densityToPressure(float &density)
{
    float density_inverse = dev_par.fluid.reference_density_inverse;
    float coef = dev_par.fluid.coefficient_rho2p;
    float gamma = dev_par.fluid.gamma;

    float pressure = (powf(density * density_inverse, gamma) - 1.0f) * coef;
    return max(pressure, dev_par.fluid.min_pressure);
}

__device__ void fluidFluidCalculate(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, float &drho_dt, glm::vec3 &du_dt, float &color)
{
    Fluid_Data fluid = dev_data.fluid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 u_i = fluid.velocity[id_i];
    glm::vec3 u_j = fluid.velocity[id_j];

    float p_i = fluid.pressure[id_i];
    float p_j = fluid.pressure[id_j];
    float cs = dev_par.fluid.sound_speed;

    float rho_i = fluid.density[id_i];
    float rho_j = fluid.density[id_j];

    glm::vec3 dx = x_j - x_i;
    float q_2 = glm::length2(dx) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    glm::vec3 dx_norm;
    if (q_2 != 0.0f)
        dx_norm = glm::normalize(dx);
    else
        dx_norm = glm::vec3(0.0f);

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = glm::dot(u_i - u_j, dx_norm);
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);

    glm::vec3 vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx_norm;

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    drho_dt += 2.0f * rho_i * glm::dot(u_i - vs, dx) * kernel_coef;

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = glm::dot(u_j - u_i, dx);
    vis /= glm::length2(dx) + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 10.0f * dev_par.fluid.viscosity * vis / rho_i * dx * kernel_coef;
}

__device__ void fluidSolidCalculate(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, float &drho_dt, glm::vec3 &du_dt, float &color)
{
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 normal = solid.normal[id_j];
    glm::vec3 solid_vel = solid.velocity[id_j];

    glm::vec3 u_i = fluid.velocity[id_i];
    glm::vec3 u_j = u_i - 2.0f * glm::dot(u_i - solid_vel, normal) * normal;

    glm::vec3 dx = x_j - x_i;

    float p_i = fluid.pressure[id_i];
    float p_j = solid.pressure[id_j];
    float cs = dev_par.fluid.sound_speed;

    float rho_i = fluid.density[id_i];
    float rho_j = rho_i;

    float q_2 = glm::length2(dx) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    glm::vec3 dx_norm;
    if (q_2 != 0.0f)
        dx_norm = glm::normalize(dx);
    else
        dx_norm = glm::vec3(0.0f);

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = glm::dot(u_i - u_j, dx_norm);
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * cs;

    glm::vec3 vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx_norm;

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    drho_dt += 2.0f * rho_i * glm::dot(u_i - vs, dx) * kernel_coef;

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = glm::dot(u_j - u_i, dx);
    vis /= glm::length2(dx) + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 10.0f * dev_par.fluid.viscosity * vis / rho_i * dx * kernel_coef;
}

__device__ void fluidVirtCalculate(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, float &drho_dt, glm::vec3 &du_dt, float &color)
{
    Fluid_Data fluid = dev_data.fluid;
    Virt_Data virt = dev_data.virt;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 u_i = fluid.velocity[id_i];
    glm::vec3 u_j = virt.velocity[id_j];

    glm::vec3 dx = x_j - x_i;
    float q_2 = glm::length2(dx) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    float rho_i = fluid.density[id_i];
    float rho_j = virt.density[id_j];

    float p_i = fluid.pressure[id_i];
    float p_j = virt.pressure[id_j];
    float cs = dev_par.fluid.sound_speed;

    glm::vec3 dx_norm;
    if (q_2 != 0.0f)
        dx_norm = glm::normalize(dx);
    else
        dx_norm = glm::vec3(0.0f);

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = glm::dot(u_i - u_j, dx_norm);
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * min(3.0f * du, cs);

    ps = max(ps, 0.0f);

    glm::vec3 vs = (rho_i * u_i + rho_j * u_j) / (rho_i + rho_j);
    vs += (p_i - p_j) / (rho_i + rho_j) / cs * dx_norm;

    du_dt += -2.0f * ps / rho_i * dx * kernel_coef;
    drho_dt += 2.0f * rho_i * glm::dot(u_i - vs, dx) * kernel_coef;

    if (!dev_par.fluid.is_fluid_viscid)
        return;
    float vis = glm::dot(u_j - u_i, dx);
    vis /= glm::length2(dx) + 0.01f * dev_par.kernel.smoothing_length_square;
    du_dt += 10.0f * dev_par.fluid.viscosity * vis / rho_i * dx * kernel_coef;
}

__global__ void findFluidPair()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Virt_Data virt = dev_data.virt;
    Cell *cell_i = dev_data.cell + fluid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    float drho_dt = 0.0f;
    glm::vec3 du_dt = glm::vec3(0.0f);
    float color = 0.0f;

    int fn = 0, vn = 0;

    glm::vec3 x_i = fluid.position[id];
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];

        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int id_j = cell_j->fluid_start_id;
        while (id_j != -1)
        {
            fn++;
            glm::vec3 x_j = fluid.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                fluidFluidCalculate(id, id_j, x_i, x_j, drho_dt, du_dt, color);
            id_j = fluid.next_id[id_j];
        }

        id_j = cell_j->solid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_j = solid.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                fluidSolidCalculate(id, id_j, x_i, x_j, drho_dt, du_dt, color);
            id_j = solid.next_id[id_j];
        }

        id_j = cell_j->virt_start_id;
        while (id_j != -1)
        {
            vn++;
            glm::vec3 x_j = virt.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                fluidVirtCalculate(id, id_j, x_i, x_j, drho_dt, du_dt, color);
            id_j = virt.next_id[id_j];
        }
    }

    fluid.drho_dt[id] = drho_dt;
    fluid.du_dt[id] = du_dt;
    fluid.color[id] = color;
}

__global__ void fluidUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.fluid.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    glm::vec3 gravity = glm::vec3(dev_par.physics.gravity[0], dev_par.physics.gravity[1], dev_par.physics.gravity[2]);
    glm::vec3 pos_min = glm::vec3(dev_par.domain.domain_min[0], dev_par.domain.domain_min[1], dev_par.domain.domain_min[2]);
    glm::vec3 pos_max = glm::vec3(dev_par.domain.domain_max[0], dev_par.domain.domain_max[1], dev_par.domain.domain_max[2]);

    float dt = dev_par.time.dt;
    if (!(fluid.drho_dt[id] == fluid.drho_dt[id]))
        fluid.drho_dt[id] = 0.0f;
    if (!(fluid.du_dt[id] == fluid.du_dt[id]))
        fluid.du_dt[id] = glm::vec3(0.0f);

    fluid.position[id] += fluid.velocity[id] * dt;
    fluid.velocity[id] += (fluid.du_dt[id] + gravity) * dt;
    fluid.density[id] += fluid.drho_dt[id] * dt;
    fluid.density[id] = min(max(fluid.density[id], dev_par.fluid.min_density), dev_par.fluid.max_density);
    fluid.pressure[id] = densityToPressure(fluid.density[id]);

    // saturate
    fluid.position[id].x = max(fluid.position[id].x, pos_min[0]);
    fluid.position[id].x = min(fluid.position[id].x, pos_max[0]);
    fluid.position[id].y = max(fluid.position[id].y, pos_min[1]);
    fluid.position[id].y = min(fluid.position[id].y, pos_max[1]);
    fluid.position[id].z = max(fluid.position[id].z, pos_min[2]);
    fluid.position[id].z = min(fluid.position[id].z, pos_max[2]);

    float speed = glm::length(fluid.velocity[id]);
    float max_speed = 0.3f * dev_par.fluid.sound_speed;
    if (speed > max_speed)
        fluid.velocity[id] *= max_speed / speed;
}

void DeviceCalculation::fluidCalculation()
{
    if (!(Sim::parameter.hasFluid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int fluid_num = Sim::parameter.fluid.number;
    int fluid_block = fluid_num / thread_num + 1;

    findFluidPair<<<fluid_block, thread_num>>>();
    cudaDeviceSynchronize();
    fluidUpdate<<<fluid_block, thread_num>>>();
    cudaDeviceSynchronize();
}