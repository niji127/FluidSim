#ifndef _PROBE_H
#define _PROBE_H

#include "sim.cuh"

#include "glm/glm.hpp"

class Probe
{
public:
    int id;
    glm::vec3 position, velocity;
    float pressure, density;
    ParticleType particle_type;
    float other_output[4];
    Probe()
    {
        id = 0;
        velocity=glm::vec3(0.0f);
        pressure = 0.0f;
        particle_type = FLUID;
        for (int i = 0; i < 4; i++)
            other_output[i] = 0.0f;
    }
    Probe(const glm::vec3 &pos, const int &probe_id) : position(pos), id(probe_id)
    {
        velocity=glm::vec3(0.0f);
        pressure = 0.0f;
        particle_type = FLUID;
        for (int i = 0; i < 4; i++)
            other_output[i] = 0.0f;
    }
    void getData();
    void outputData();
};

#endif