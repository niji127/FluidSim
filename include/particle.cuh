#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "global.h"
#include "types.h"
#include "parameter.h"
#include "host_data.cuh"

#include "glm/glm.hpp"
#include "thrust/device_vector.h"
#include "thrust/sort.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/tuple.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <string>
#include <sstream>

class Fluid_List
{
public:
    thrust::device_vector<glm::vec3> position, velocity;
    thrust::device_vector<float> density, pressure, color;
    thrust::device_vector<int> cell_id, next_id, next_id_previous;

    thrust::device_vector<glm::vec3> du_dt;
    thrust::device_vector<float> drho_dt;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += position.size() * sizeof(glm::vec3) + sizeof(position);
        total_size += velocity.size() * sizeof(glm::vec3) + sizeof(velocity);
        total_size += density.size() * sizeof(float) + sizeof(density);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += color.size() * sizeof(float) + sizeof(color);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);
        total_size += next_id.size() * sizeof(int) + sizeof(next_id);
        total_size += next_id_previous.size() * sizeof(int) + sizeof(next_id_previous);

        total_size += du_dt.size() * sizeof(glm::vec3) + sizeof(du_dt);
        total_size += drho_dt.size() * sizeof(float) + sizeof(drho_dt);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Fluid_Host_List &host_data)
    {
        position = host_data.position;
        velocity = host_data.velocity;
        density = host_data.density;
        pressure = host_data.pressure;
        color = host_data.color;
        cell_id = host_data.cell_id;
        next_id = host_data.next_id;
        next_id_previous = host_data.next_id_previous;

        du_dt = host_data.du_dt;
        drho_dt = host_data.drho_dt;
    }

    void initiate(const Parameter &parameter)
    {
        int fluid_num = parameter.fluid.number;

        cell_id_list.reserve(fluid_num);
        cell_counts.reserve(fluid_num);
        cell_indices.reserve(fluid_num);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(position.begin(), velocity.begin(), density.begin(), pressure.begin(), color.begin(),next_id_previous.begin()));
    }

    void sort();
};

class Solid_List
{
public:
    thrust::device_vector<glm::vec3> position, reference_position, velocity;
    thrust::device_vector<SolidType> solid_type;
    thrust::device_vector<float> pressure;
    thrust::device_vector<int> cell_id, next_id, next_id_previous;

    thrust::device_vector<glm::mat3> deformation, stress_1st_piola, correction;
    thrust::device_vector<glm::vec3> du_dt;
    thrust::device_vector<glm::vec3> normal;
    thrust::device_vector<glm::vec3> couple_dudt;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += reference_position.size() * sizeof(glm::vec3) + sizeof(reference_position);
        total_size += position.size() * sizeof(glm::vec3) + sizeof(position);
        total_size += velocity.size() * sizeof(glm::vec3) + sizeof(velocity);
        total_size += solid_type.size() * sizeof(SolidType) + sizeof(solid_type);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);
        total_size += next_id.size() * sizeof(int) + sizeof(next_id);
        total_size += next_id_previous.size() * sizeof(int) + sizeof(next_id_previous);

        total_size += deformation.size() * sizeof(glm::mat3) + sizeof(deformation);
        total_size += stress_1st_piola.size() * sizeof(glm::mat3) + sizeof(stress_1st_piola);
        total_size += correction.size() * sizeof(glm::mat3) + sizeof(correction);
        total_size += du_dt.size() * sizeof(glm::vec3) + sizeof(du_dt);
        total_size += normal.size() * sizeof(glm::vec3) + sizeof(normal);
        total_size += couple_dudt.size() * sizeof(glm::vec3) + sizeof(couple_dudt);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Solid_Host_List &host_data)
    {
        position = host_data.position;
        reference_position = host_data.reference_position;
        velocity = host_data.velocity;
        solid_type = host_data.solid_type;
        pressure = host_data.pressure;
        cell_id = host_data.cell_id;
        next_id = host_data.next_id;
        next_id_previous = host_data.next_id_previous;

        deformation = host_data.deformation;
        stress_1st_piola = host_data.stress_1st_piola;
        correction = host_data.correction;
        du_dt = host_data.du_dt;
        normal = host_data.normal;
        couple_dudt = host_data.couple_dudt;
    }

    void initiate(const Parameter &parameter)
    {
        cell_id_list.reserve(parameter.solid.number);
        cell_counts.reserve(parameter.solid.number);
        cell_indices.reserve(parameter.solid.number);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(reference_position.begin(), position.begin(), velocity.begin(), solid_type.begin(), du_dt.begin(),next_id_previous.begin()));
    }

    void sort();
};

class Virt_List
{
public:
    thrust::device_vector<glm::vec3> position, normal, velocity;
    thrust::device_vector<float> pressure, density;
    thrust::device_vector<VirtType> virt_type;
    thrust::device_vector<int> cell_id, next_id, next_id_previous;

    thrust::device_vector<int> cell_id_list;
    thrust::device_vector<int> cell_counts;
    thrust::device_vector<int> cell_indices;
    int cell_num;

    size_t getMemoryUsed()
    {
        size_t total_size = 0;
        total_size += position.size() * sizeof(glm::vec3) + sizeof(position);
        total_size += normal.size() * sizeof(glm::vec3) + sizeof(normal);
        total_size += velocity.size() * sizeof(glm::vec3) + sizeof(velocity);
        total_size += pressure.size() * sizeof(float) + sizeof(pressure);
        total_size += density.size() * sizeof(float) + sizeof(density);
        total_size += virt_type.size() * sizeof(VirtType) + sizeof(virt_type);
        total_size += cell_id.size() * sizeof(int) + sizeof(cell_id);
        total_size += next_id.size() * sizeof(int) + sizeof(next_id);
        total_size += next_id_previous.size() * sizeof(int) + sizeof(next_id_previous);

        total_size += cell_id_list.size() * sizeof(int) + sizeof(cell_id_list);
        total_size += cell_counts.size() * sizeof(int) + sizeof(cell_counts);
        total_size += cell_indices.size() * sizeof(int) + sizeof(cell_indices);
        total_size += sizeof(cell_num);

        return total_size;
    }

    void copyData(const Virt_Host_List &host_data)
    {
        position = host_data.position;
        normal = host_data.normal;
        velocity = host_data.velocity;
        pressure = host_data.pressure;
        density = host_data.density;
        virt_type = host_data.virt_type;
        cell_id = host_data.cell_id;
        next_id = host_data.next_id;
        next_id_previous = host_data.next_id_previous;
    }

    void initiate(const Parameter &parameter)
    {
        cell_id_list.reserve(parameter.virt.number);
        cell_counts.reserve(parameter.virt.number);
        cell_indices.reserve(parameter.virt.number);
        cell_num = 0;
    }

    auto getZipProperties()
    {
        return thrust::make_zip_iterator(thrust::make_tuple(position.begin(), virt_type.begin(),next_id_previous.begin()));
    }

    void sort();
};
#endif