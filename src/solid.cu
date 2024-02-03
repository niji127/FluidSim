#include "calculation.cuh"
#include "memory_copy.cuh"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>

extern __device__ float densityToPressure(float &density);
extern __device__ float pressureToDensity(float &pressure);

__device__ void solidSolidDeformation(int id_i, int id_j, const glm::vec3 &x_ref_i, const glm::vec3 &x_ref_j, glm::mat3 &deformation, glm::mat3 &correction)
{
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 dx = solid.position[id_j] - solid.position[id_i];

    glm::vec3 dx_ref = x_ref_j - x_ref_i;
    float q_2 = glm::length2(dx_ref) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    deformation += glm::outerProduct(dx, dx_ref) * kernel_coef;
    correction += glm::outerProduct(dx_ref, dx_ref) * kernel_coef;
}

__device__ glm::mat3 solidStressCalculation(const glm::mat3 &deformation)
{
    float lambda = dev_par.solid.lambda;
    float nu = dev_par.solid.nu;

    glm::mat3 green_strain = 0.5f * (glm::transpose(deformation) * deformation - glm::mat3(1.0f));
    float trace = green_strain[0][0] + green_strain[1][1] + green_strain[2][2];
    glm::mat3 stress_2nd_piola = 2.0f * nu * green_strain + lambda * trace * glm::mat3(1.0f);
    return deformation * stress_2nd_piola;
}

__global__ void solidDeformationCalculate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;

    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    glm::mat3 correction = glm::mat3(0.0f);
    glm::mat3 deformation = glm::mat3(0.0f);

    glm::vec3 x_ref_i = solid.reference_position[id];
    // deformation && correction
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int id_j = cell_j->solid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_ref_j = solid.reference_position[id_j];
            float dis_ref_2 = glm::length2(x_ref_j - x_ref_i);
            if (dis_ref_2 < impact_dis_2)
                solidSolidDeformation(id, id_j, x_ref_i, x_ref_j, deformation, correction);
            id_j = solid.next_id[id_j];
        }
    }

    if (glm::determinant(correction) != 0.0f)
        correction = glm::inverse(correction);
    else
        correction = glm::mat3(1.0f);
    deformation = deformation * correction;
    solid.deformation[id] = deformation;
    solid.stress_1st_piola[id] = solidStressCalculation(deformation);
    solid.correction[id] = correction;
}

__device__ void solidFluidCouple(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, glm::vec3 &du_dt)
{
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 normal = solid.normal[id_i];
    glm::vec3 solid_vel = solid.velocity[id_i];

    glm::vec3 u_j = fluid.velocity[id_j];
    glm::vec3 u_i = u_j - 2.0f * glm::dot(u_j - solid_vel, normal) * normal;

    glm::vec3 dx = x_j - x_i;
    float q_2 = glm::length2(dx) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    float p_i = solid.pressure[id_i];
    float p_j = fluid.pressure[id_j];
    float cs = dev_par.fluid.sound_speed;

    float rho_j = fluid.density[id_j];
    float rho_i = rho_j;

    glm::vec3 dx_norm = glm::vec3(0.0f);
    if (glm::length(dx) != 0.0f)
        dx_norm = glm::normalize(dx);

    float ps = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
    float du = glm::dot(u_i - u_j, dx_norm);
    if (du >= 0.0)
        ps += rho_i * rho_j * du / (rho_i + rho_j) * cs;

    float rho_ref = dev_par.solid.reference_density;
    du_dt += -2.0f * ps / rho_ref * dx * kernel_coef;
}

__device__ void solidSolidForce(int id_i, int id_j, glm::vec3 &x_ref_i, glm::vec3 &x_ref_j, glm::vec3 &du_dt)
{
    Solid_Data solid = dev_data.solid;

    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;
    float rho_ref = dev_par.solid.reference_density;

    glm::vec3 dx_ref = x_ref_j - x_ref_i;
    float q_2 = glm::length2(dx_ref) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    glm::mat3 deformation_i = solid.deformation[id_i];
    glm::mat3 stress_i = solid.stress_1st_piola[id_i];
    glm::mat3 stress_j = solid.stress_1st_piola[id_j];
    glm::mat3 correction_i = solid.correction[id_i];
    glm::mat3 correction_j = solid.correction[id_j];

    du_dt += (stress_j * correction_j + stress_i * correction_i) * dx_ref * kernel_coef / rho_ref;

    glm::vec3 dx = solid.position[id_j] - solid.position[id_i];
    glm::vec3 du = solid.velocity[id_j] - solid.velocity[id_i];

    float du_dot_dx = glm::dot(du, dx);

    float alpha = dev_par.solid.artificial_viscocity[0];
    float beta = dev_par.solid.artificial_viscocity[1];
    float eps = 0.01f * dev_par.kernel.smoothing_length_square;
    float cs = dev_par.solid.sound_speed;

    float dis_square = glm::length2(dx);
    float nu = dev_par.kernel.smoothing_length * du_dot_dx / (dis_square + eps);

    float viscosity = alpha * nu * cs;
    viscosity -= beta * nu * nu;

    float det = glm::determinant(deformation_i);
    du_dt += det * glm::transpose(glm::inverse(deformation_i)) * viscosity * correction_i * dx_ref * kernel_coef;
}

__global__ void solidForceCalculate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;

    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;

    glm::vec3 couple_du_dt = glm::vec3(0.0f);
    glm::vec3 du_dt = glm::vec3(0.0f);

    glm::vec3 x_i = solid.position[id];
    glm::vec3 x_ref_i = solid.reference_position[id];

    // du_dt
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int id_j = cell_j->fluid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_j = fluid.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                solidFluidCouple(id, id_j, x_i, x_j, couple_du_dt);
            id_j = fluid.next_id[id_j];
        }

        id_j = cell_j->solid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_ref_j = solid.reference_position[id_j];
            float dis_ref_2 = glm::length2(x_ref_j - x_ref_i);
            if (dis_ref_2 < impact_dis_2)
                solidSolidForce(id, id_j, x_ref_i, x_ref_j, du_dt);
            id_j = solid.next_id[id_j];
        }
    }
    solid.couple_dudt[id] = couple_du_dt;
    solid.du_dt[id] = du_dt + couple_du_dt;
}

__global__ void solidUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;
    Solid_Data solid = dev_data.solid;
    glm::vec3 gravity = glm::vec3(dev_par.physics.gravity[0], dev_par.physics.gravity[1], dev_par.physics.gravity[2]);
    glm::vec3 pos_min = glm::vec3(dev_par.domain.domain_min[0], dev_par.domain.domain_min[1], dev_par.domain.domain_min[2]);
    glm::vec3 pos_max = glm::vec3(dev_par.domain.domain_max[0], dev_par.domain.domain_max[1], dev_par.domain.domain_max[2]);

    float dt = dev_par.time.dt / (float)dev_par.time.solid_sub_step;
    switch (solid.solid_type[id])
    {
    case FIXED_SOLID:
        solid.velocity[id] = glm::vec3(0.0f);
        solid.position[id] += solid.velocity[id] * dt;
        solid.du_dt[id] = glm::vec3(0.0f);
        break;
    case UNCONSTRAINED_SOLID:
        solid.position[id] += solid.velocity[id] * dt;
        solid.velocity[id] += (solid.du_dt[id] + gravity) * dt;
        // saturate
        solid.position[id].x = max(solid.position[id].x, pos_min[0]);
        solid.position[id].x = min(solid.position[id].x, pos_max[0]);
        solid.position[id].y = max(solid.position[id].y, pos_min[1]);
        solid.position[id].y = min(solid.position[id].y, pos_max[1]);
        solid.position[id].z = max(solid.position[id].z, pos_min[2]);
        solid.position[id].z = min(solid.position[id].z, pos_max[2]);
        break;
    }
}

void DeviceCalculation::solidCalculation(Sim *sim)
{
    if (!(Sim::parameter.hasSolid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;

    for (int i = 0; i < Sim::parameter.time.solid_sub_step; i++)
    {
        SimTime *time = Sim::parameter.getTime();
        solidDeformationCalculate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
        solidForceCalculate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
        solidUpdate<<<solid_block, thread_num>>>();
        cudaDeviceSynchronize();
    }
}

__device__ void solidNormalCalculate(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, glm::vec3 &normal)
{
    float coef = dev_par.kernel.kernel_differential;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 dx = x_j - x_i;

    float q_2 = glm::length2(dx) * h_inv_2;
    float kernel_coef = expf(-q_2) * coef;

    normal -= dx * kernel_coef;
}

__device__ void solidPressureCalculate(int id_i, int id_j, glm::vec3 &x_i, glm::vec3 &x_j, float &pressure, float &sum_correction)
{
    Fluid_Data fluid = dev_data.fluid;

    float coef1 = dev_par.kernel.kernel_coefficient_1;
    float coef2 = dev_par.kernel.kernel_coefficient_2;
    float h_inv_2 = dev_par.kernel.smoothing_length_inverse_square;

    glm::vec3 dx = x_j - x_i;

    float q_2 = glm::length2(dx) * h_inv_2;
    float pre_j = fluid.pressure[id_j];
    float kernel_coef = (expf(-q_2) - coef1) * coef2;

    sum_correction += kernel_coef;
    pressure += pre_j * kernel_coef;
}

__global__ void solidCoupleUpdate()
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= dev_par.solid.number)
        return;
    Fluid_Data fluid = dev_data.fluid;
    Solid_Data solid = dev_data.solid;
    Cell *cell_i = dev_data.cell + solid.cell_id[id];

    float impact_dis_2 = dev_par.kernel.impact_distance_square;
    glm::vec3 normal = glm::vec3(0.0f);
    float pressure(0.0f), sum_correction(0.0f);

    glm::vec3 x_i = solid.position[id];
    // normal
    for (int n_id = 0; n_id < 27; n_id++)
    {
        int neighbor_id = cell_i->neighbor_cell_id[n_id];
        if (neighbor_id == -1)
            continue;
        Cell *cell_j = dev_data.cell + neighbor_id;

        int id_j = cell_j->fluid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_j = fluid.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                solidPressureCalculate(id, id_j, x_i, x_j, pressure, sum_correction);
            id_j = fluid.next_id[id_j];
        }

        id_j = cell_j->solid_start_id;
        while (id_j != -1)
        {
            glm::vec3 x_j = solid.position[id_j];
            float dis_2 = glm::length2(x_j - x_i);
            if (dis_2 < impact_dis_2)
                solidNormalCalculate(id, id_j, x_i, x_j, normal);
            id_j = solid.next_id[id_j];
        }
    }

    if (glm::length2(normal) < EPS_FOR_SUM)
        solid.normal[id] = glm::vec3(0.0f);
    else
        solid.normal[id] = glm::normalize(normal);
    if (sum_correction > EPS_FLOAT)
        solid.pressure[id] = max(pressure / sum_correction, dev_par.fluid.min_pressure);
    else
        solid.pressure[id] = 0.0f;
}

void DeviceCalculation::coupleCalculation()
{
    if (!(Sim::parameter.hasFluid()))
        return;
    if (!(Sim::parameter.hasSolid()))
        return;

    int thread_num = Sim::parameter.kernel.thread_num;
    int solid_num = Sim::parameter.solid.number;
    int solid_block = solid_num / thread_num + 1;

    solidCoupleUpdate<<<solid_block, thread_num>>>();
    cudaDeviceSynchronize();
}