#include "properties.h"
#include <iostream>
#include <fstream>
#include <string>
#include <windows.h>

bool FluidProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\fluid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read fluid particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number;
    if (number == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number)
    {
        std::cout << "mismatch fluid particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool FluidProperties::readParameter(const Config &config)
{
    is_fluid_viscid = config.Read("is_viscid", false);
    min_pressure = config.Read("min_pressure", -1.013e5f);

    sound_speed = config.Read("fluid_sound_speed", 1500.0f);
    reference_density = config.Read("fluid_density", 1000.0f);
    gamma = config.Read("gamma", 7.0f);
    viscosity = config.Read("viscosity", 1.01e-3f); // 20 'C

    return false;
}

void FluidProperties::initiate()
{
    float variation_percent = 0.1f;
    coefficient_p2rho = gamma * pow(sound_speed, -2.0f) / reference_density;
    coefficient_rho2p = 1.0f / coefficient_p2rho;
    gamma_inv = 1.0f / gamma;
    reference_density_inverse = 1.0f / reference_density;

    min_density = (1.0f - variation_percent) * reference_density;
    max_density = (1.0f + variation_percent) * reference_density;
}

void SolidProperties::initiate()
{
    lambda = youngs_modulus * poisson_ratio / (1.0f - 2.0f * poisson_ratio) / (1.0f + poisson_ratio);
    nu = youngs_modulus / (1.0f + poisson_ratio) / 2.0f;

    bulk_modulus = youngs_modulus / (3.0f * (1.0f - 2.0f * poisson_ratio));
    sound_speed = sqrtf(bulk_modulus / reference_density);
}

bool SolidProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\solid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read solid particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number;
    if (number == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number)
    {
        std::cout << "mismatch solid particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool SolidProperties::readParameter(const Config &config)
{
    reference_density = config.Read("solid_density", 2700.0f);
    youngs_modulus = config.Read("youngs_modulus", 67.5e9f);
    poisson_ratio = config.Read("poisson_ratio", 0.34f);
    artificial_viscocity[0] = config.Read("arti_vis_alpha", 2.5f);
    artificial_viscocity[1] = config.Read("arti_vis_beta", 2.5f);

    return false;
}

void VirtualProperties::initiate() {}
bool VirtualProperties::readNumber(const std::string &path)
{
    std::ifstream infile;
    infile.open(path + "\\virt_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read virtual particles parameter" << std::endl;
        return true;
    }
    std::string str;
    std::getline(infile, str);
    std::istringstream string_to_num(str);
    string_to_num >> number;
    if (number == 0)
        return false;
    int lines_check(0);
    while (std::getline(infile, str))
        lines_check++;
    if (lines_check != number)
    {
        std::cout << "mismatch virtual particle numbers" << std::endl;
        return true;
    }
    infile.close();
    return false;
}

bool VirtualProperties::readParameter(const Config &config)
{
    return false;
}
