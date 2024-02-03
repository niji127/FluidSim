#ifndef _DATA_POINTER_CUH
#define _DATA_POINTER_CUH

#include "particle.cuh"
#include "cell.cuh"

class Fluid_Data
{
public:
    glm::vec3 *position, *velocity;
    float *density, *pressure, *color;
    int *cell_id, *next_id, *next_id_previous;
    glm::vec3 *du_dt;
    float *drho_dt;
};

class Solid_Data
{
public:
    glm::vec3 *position, *reference_position, *velocity;
    SolidType *solid_type;
    float *pressure;
    int *cell_id, *next_id, *next_id_previous;

    glm::mat3 *deformation, *stress_1st_piola, *correction;
    glm::vec3 *du_dt;
    glm::vec3 *couple_dudt;
    glm::vec3 *normal;
};

class Virt_Data
{
public:
    glm::vec3 *position, *normal, *velocity;
    float *pressure, *density;
    VirtType *virt_type;
    int *cell_id, *next_id, *next_id_previous;
};

class Data_Pointer
{
public:
    Fluid_Data fluid;
    Solid_Data solid;
    Virt_Data virt;
    Cell *cell;
    void fluidInitiate(Fluid_List &fluid_list)
    {
        fluid.position = thrust::raw_pointer_cast(fluid_list.position.data());
        fluid.velocity = thrust::raw_pointer_cast(fluid_list.velocity.data());
        fluid.density = thrust::raw_pointer_cast(fluid_list.density.data());
        fluid.pressure = thrust::raw_pointer_cast(fluid_list.pressure.data());
        fluid.color = thrust::raw_pointer_cast(fluid_list.color.data());
        fluid.cell_id = thrust::raw_pointer_cast(fluid_list.cell_id.data());
        fluid.next_id = thrust::raw_pointer_cast(fluid_list.next_id.data());
        fluid.next_id_previous = thrust::raw_pointer_cast(fluid_list.next_id_previous.data());

        fluid.du_dt = thrust::raw_pointer_cast(fluid_list.du_dt.data());
        fluid.drho_dt = thrust::raw_pointer_cast(fluid_list.drho_dt.data());
    }

    void solidInitiate(Solid_List &solid_list)
    {
        solid.reference_position = thrust::raw_pointer_cast(solid_list.reference_position.data());
        solid.position = thrust::raw_pointer_cast(solid_list.position.data());
        solid.velocity = thrust::raw_pointer_cast(solid_list.velocity.data());
        solid.solid_type = thrust::raw_pointer_cast(solid_list.solid_type.data());
        solid.pressure = thrust::raw_pointer_cast(solid_list.pressure.data());
        solid.cell_id = thrust::raw_pointer_cast(solid_list.cell_id.data());
        solid.next_id = thrust::raw_pointer_cast(solid_list.next_id.data());
        solid.next_id_previous = thrust::raw_pointer_cast(solid_list.next_id_previous.data());

        solid.deformation = thrust::raw_pointer_cast(solid_list.deformation.data());
        solid.stress_1st_piola = thrust::raw_pointer_cast(solid_list.stress_1st_piola.data());
        solid.correction = thrust::raw_pointer_cast(solid_list.correction.data());
        solid.du_dt = thrust::raw_pointer_cast(solid_list.du_dt.data());
        solid.couple_dudt = thrust::raw_pointer_cast(solid_list.couple_dudt.data());
        solid.normal = thrust::raw_pointer_cast(solid_list.normal.data());
    }

    void virtInitiate(Virt_List &virt_list)
    {
        virt.position = thrust::raw_pointer_cast(virt_list.position.data());
        virt.normal = thrust::raw_pointer_cast(virt_list.normal.data());
        virt.velocity = thrust::raw_pointer_cast(virt_list.velocity.data());
        virt.pressure = thrust::raw_pointer_cast(virt_list.pressure.data());
        virt.density = thrust::raw_pointer_cast(virt_list.density.data());
        virt.virt_type = thrust::raw_pointer_cast(virt_list.virt_type.data());
        virt.cell_id = thrust::raw_pointer_cast(virt_list.cell_id.data());
        virt.next_id = thrust::raw_pointer_cast(virt_list.next_id.data());
        virt.next_id_previous = thrust::raw_pointer_cast(virt_list.next_id_previous.data());
    }

    void cellInitiate(thrust::device_vector<Cell> &cell_list)
    {
        cell = thrust::raw_pointer_cast(cell_list.data());
    }
};

#endif