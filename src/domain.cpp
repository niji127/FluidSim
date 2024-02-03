#include "domain.h"
#include "global.h"

bool Domain::readNumber(const std::string &path)
{
    return false;
}

bool Domain::readParameter(const Config &config)
{
    domain_min[0] = config.Read("xmin", 0.0f);
    domain_min[1] = config.Read("ymin", 0.0f);
    domain_min[2] = config.Read("zmin", 0.0f);

    domain_max[0] = config.Read("xmax", 0.0f);
    domain_max[1] = config.Read("ymax", 0.0f);
    domain_max[2] = config.Read("zmax", 0.0f);

    if (domain_max[0] - domain_min[0] == 0.0f || domain_max[1] - domain_min[1] == 0.0f || domain_max[2] - domain_min[2] == 0.0f)
    {
        std::cout << "failed to read domain information" << std::endl;
        return true;
    }

    return false;
}

void Domain::initiate(const Kernel &kernel)
{
    float boudary_thick = (float)VIRT_LAYER * kernel.particle_diameter;
    for (int i = 0; i < 3; i++)
    {
        cell_min[i] = domain_min[i] - boudary_thick;
        cell_max[i] = domain_max[i] + boudary_thick;
    }

    float smoothing_length = kernel.smoothing_length;
    float cell_size_max = kernel.impact_length_hsml_ratio * smoothing_length;

    for (int i = 0; i < 3; i++)
    {
        cell_number[i] = (int)((cell_max[i] - cell_min[i]) / cell_size_max);
        interval[i] = (cell_max[i] - cell_min[i]) / (float)cell_number[i];
        interval_inv[i] = 1.0f / interval[i];
    }
    cell_number_total = cell_number[0] * cell_number[1] * cell_number[2];
}