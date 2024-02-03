#ifndef _HOST_DATA_CUH
#define _HOST_DATA_CUH

#include "parameter.h"
#include "types.h"

#include "glm/glm.hpp"
#include "thrust/host_vector.h"

class Fluid_Host_List
{
public:
    thrust::host_vector<glm::vec3> position, velocity;
    thrust::host_vector<float> density, pressure, color;
    thrust::host_vector<int> cell_id, next_id, next_id_previous;

    thrust::host_vector<glm::vec3> du_dt;
    thrust::host_vector<float> drho_dt;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);
        glm::vec3 temp_pos, temp_velo;
        float temp_dens, temp_pres;
        string_to_num >> temp_pos.x >> temp_pos.y >> temp_pos.z >> temp_velo.x >> temp_velo.y >> temp_velo.z >> temp_dens >> temp_pres;

        position.push_back(temp_pos);
        velocity.push_back(temp_velo);
        density.push_back(temp_dens);
        pressure.push_back(temp_pres);
        color.push_back(0.0f);
        cell_id.push_back(-1);
        next_id.push_back(-1);
        next_id_previous.push_back(-1);

        du_dt.push_back(glm::vec3(0.0f));
        drho_dt.push_back(0.0f);
    }

    void clear()
    {
        position.clear();
        velocity.clear();
        density.clear();
        pressure.clear();
        color.clear();
        cell_id.clear();
        next_id.clear();
        next_id_previous.clear();
        du_dt.clear();
        drho_dt.clear();

        position.shrink_to_fit();
        velocity.shrink_to_fit();
        density.shrink_to_fit();
        pressure.shrink_to_fit();
        color.shrink_to_fit();
        cell_id.shrink_to_fit();
        next_id.shrink_to_fit();
        next_id_previous.shrink_to_fit();
        du_dt.shrink_to_fit();
        drho_dt.shrink_to_fit();
    }
};

class Solid_Host_List
{
public:
    thrust::host_vector<glm::vec3> position, reference_position, velocity;
    thrust::host_vector<SolidType> solid_type;
    thrust::host_vector<float> pressure;
    thrust::host_vector<int> cell_id, next_id, next_id_previous;

    thrust::host_vector<glm::mat3> deformation, stress_1st_piola, correction;
    thrust::host_vector<glm::vec3> du_dt;
    thrust::host_vector<glm::vec3> normal;
    thrust::host_vector<glm::vec3> couple_dudt;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);
        int solid_type_id;
        glm::vec3 temp_pos, temp_ref_pos, temp_velo;
        string_to_num >> temp_pos.x >> temp_pos.y >> temp_pos.z >> temp_velo.x >> temp_velo.y >> temp_velo.z >> temp_ref_pos.x >> temp_ref_pos.y >> temp_ref_pos.z >> solid_type_id;

        reference_position.push_back(temp_ref_pos);
        position.push_back(temp_pos);
        velocity.push_back(temp_velo);
        solid_type.push_back(SolidType(solid_type_id));
        pressure.push_back(0.0f);
        cell_id.push_back(-1);
        next_id.push_back(-1);
        next_id_previous.push_back(-1);

        deformation.push_back(glm::mat3(1.0f));
        stress_1st_piola.push_back(glm::mat3(0.0f));
        correction.push_back(glm::mat3(1.0f));
        du_dt.push_back(glm::vec3(0.0f));
        couple_dudt.push_back(glm::vec3(0.0f));
        normal.push_back(glm::vec3(0.0f));
    }

    void clear()
    {
        position.clear();
        reference_position.clear();
        velocity.clear();
        solid_type.clear();
        pressure.clear();
        cell_id.clear();
        next_id.clear();
        next_id_previous.clear();

        deformation.clear();
        stress_1st_piola.clear();
        correction.clear();
        du_dt.clear();
        normal.clear();
        couple_dudt.clear();

        position.shrink_to_fit();
        reference_position.shrink_to_fit();
        velocity.shrink_to_fit();
        solid_type.shrink_to_fit();
        pressure.shrink_to_fit();
        cell_id.shrink_to_fit();
        next_id.shrink_to_fit();
        next_id_previous.shrink_to_fit();

        deformation.shrink_to_fit();
        stress_1st_piola.shrink_to_fit();
        correction.shrink_to_fit();
        du_dt.shrink_to_fit();
        normal.shrink_to_fit();
        couple_dudt.shrink_to_fit();
    }
};

class Virt_Host_List
{
public:
    thrust::host_vector<glm::vec3> position, normal, velocity;
    thrust::host_vector<float> pressure, density;
    thrust::host_vector<VirtType> virt_type;
    thrust::host_vector<int> cell_id, next_id, next_id_previous;

    void addList(const std::string &str, const Parameter &parameter)
    {
        std::istringstream string_to_num(str);
        int virt_type_id;
        glm::vec3 temp_pos;
        string_to_num >> temp_pos.x >> temp_pos.y >> temp_pos.z >> virt_type_id;

        position.push_back(temp_pos);
        normal.push_back(glm::vec3(0.0f));
        velocity.push_back(glm::vec3(0.0f));
        pressure.push_back(0.0f);
        density.push_back(parameter.fluid.reference_density);
        virt_type.push_back(VirtType(virt_type_id));
        cell_id.push_back(-1);
        next_id.push_back(-1);
        next_id_previous.push_back(-1);
    }

    void clear()
    {
        position.clear();
        normal.clear();
        velocity.clear();
        pressure.clear();
        density.clear();
        virt_type.clear();
        cell_id.clear();
        next_id.clear();
        next_id_previous.clear();

        position.shrink_to_fit();
        normal.shrink_to_fit();
        velocity.shrink_to_fit();
        pressure.shrink_to_fit();
        density.shrink_to_fit();
        virt_type.shrink_to_fit();
        cell_id.shrink_to_fit();
        next_id.shrink_to_fit();
        next_id_previous.shrink_to_fit();
    }
};

#endif