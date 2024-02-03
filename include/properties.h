#ifndef _PROPERTIES_H
#define _PROPERTIES_H

#include "config.h"
#include "types.h"

class FluidProperties
{
public:
    int number;
    float sound_speed, reference_density;
    float reference_density_inverse;
    float gamma, viscosity; // dynamic viscosity
    float gamma_inv;
    float coefficient_p2rho;
    float coefficient_rho2p;
    float min_pressure;
    float min_density, max_density;
    bool is_fluid_viscid;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

class SolidProperties
{
public:
    int number;
    float reference_density, youngs_modulus, poisson_ratio;
    float lambda, nu, bulk_modulus;
    float artificial_viscocity[2], sound_speed;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

class VirtualProperties
{
public:
    int number;

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

#endif