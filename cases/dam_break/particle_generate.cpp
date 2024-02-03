#include "glm/glm.hpp"
#include "types.h"
#include "global.h"

#include <vector>
#include <fstream>
#include <iostream>

class Fluid_Particle
{
public:
    glm::vec3 position, velocity;
    float density, pressure;
};

class Solid_Particle
{
public:
    glm::vec3 velocity, position;
    glm::vec3 reference_position;
    SolidType solid_type;
};

class Virtual_Particle
{
public:
    glm::vec3 position;
    VirtType virt_type;
};

int main()
{
    float dx = 1e-2;
    float xmin(0.0f), xmax(5.366f), ymin(0.0f), ymax(2.5f), zmin(0.0f), zmax(2.0f);
    glm::vec2 fluid_boundary(2.0f, 1.0f);
    glm::vec3 solid_min = glm::vec3(4.0f, -1000.0f, 0.75f);
    glm::vec3 solid_max = glm::vec3(solid_min.x + 0.05f, 1.0f, solid_min.z + 0.5f);

    std::vector<Fluid_Particle> fluid_list;
    std::vector<Solid_Particle> solid_list;
    std::vector<Virtual_Particle> virt_list;

    int bound_layer = VIRT_LAYER;
    int nx = (int)((xmax - xmin) / dx);
    int ny = (int)((ymax - ymin) / dx);
    int nz = (int)((zmax - zmin) / dx);

    float fluid_density = 1.0;

    for (int i = -bound_layer; i < nx + bound_layer; i++)
    {
        for (int j = -bound_layer; j < ny + bound_layer; j++)
        {
            for (int k = -bound_layer; k < nz + bound_layer; k++)
            {
                glm::vec3 pos = glm::vec3(((float)i + 0.5f) * dx + xmin, ((float)j + 0.5f) * dx + ymin, ((float)k + 0.5f) * dx + zmin);

                // virtual particles
                if (pos.x < solid_max.x && pos.x > solid_min.x && pos.y < solid_max.y && pos.y > solid_min.y && pos.z < solid_max.z && pos.z > solid_min.z)
                {
                    Solid_Particle solid;
                    solid.position = pos;
                    solid.reference_position = pos;
                    solid.velocity = glm::vec3(0.0f);
                    if (pos.y < ymin)
                        solid.solid_type = FIXED_SOLID;
                    else
                        solid.solid_type = UNCONSTRAINED_SOLID;
                    solid_list.push_back(solid);
                }
                else if (pos.x < xmin || pos.x > xmax || pos.y < ymin || pos.y > ymax || pos.z < zmin || pos.z > zmax)
                {

                    Virtual_Particle virt;
                    virt.position = pos;
                    virt.virt_type = FIXED_VIRT;
                    virt_list.push_back(virt);
                }
                // liquid particles
                else if (pos.x < fluid_boundary.x && pos.y < fluid_boundary.y)
                {
                    Fluid_Particle fluid;
                    fluid.position = pos;
                    fluid.velocity = glm::vec3(0.0f);
                    fluid.density = fluid_density;
                    fluid.pressure = 0.0f;
                    fluid_list.push_back(fluid);
                }
            }
        }
    }

    std::fstream fluid_file;
    fluid_file.open("data\\fluid_particle.dat", std::ios::out);
    fluid_file << fluid_list.size() << std::endl;
    for (auto fp : fluid_list)
    {
        fluid_file << fp.position.x << " ";
        fluid_file << fp.position.y << " ";
        fluid_file << fp.position.z << " ";
        fluid_file << fp.velocity.x << " ";
        fluid_file << fp.velocity.y << " ";
        fluid_file << fp.velocity.z << " ";
        fluid_file << fp.density << " ";
        fluid_file << fp.pressure << std::endl;
    }
    fluid_file.close();

    std::fstream virt_file;
    virt_file.open("data\\virt_particle.dat", std::ios::out);
    virt_file << virt_list.size() << std::endl;
    for (auto vp : virt_list)
    {
        virt_file << vp.position.x << " ";
        virt_file << vp.position.y << " ";
        virt_file << vp.position.z << " ";
        virt_file << vp.virt_type << std::endl;
    }
    virt_file.close();

    std::fstream solid_file;
    solid_file.open("data\\solid_particle.dat", std::ios::out);
    solid_file << solid_list.size() << std::endl;
    for (auto sp : solid_list)
    {
        solid_file << sp.position.x << " ";
        solid_file << sp.position.y << " ";
        solid_file << sp.position.z << " ";
        solid_file << sp.velocity.x << " ";
        solid_file << sp.velocity.y << " ";
        solid_file << sp.velocity.z << " ";
        solid_file << sp.reference_position.x << " ";
        solid_file << sp.reference_position.y << " ";
        solid_file << sp.reference_position.z << " ";
        solid_file << sp.solid_type << std::endl;
    }
    solid_file.close();

    return 0;
}