#ifndef _DOMAIN_H
#define _DOMAIN_H

#include "config.h"
#include "kernel.h"

class Domain
{
public:
    float cell_min[3], cell_max[3];
    float domain_min[3], domain_max[3];
    int cell_number[3];
    float interval[3];
    float interval_inv[3];
    int cell_number_total;

    bool readParameter(const Config &config);
    void initiate(const Kernel &Kernel);
    bool readNumber(const std::string &path);
};

#endif