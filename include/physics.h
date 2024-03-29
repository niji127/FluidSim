#ifndef _PHYSICS_H
#define _PHYSICS_H

#include "config.h"
#include "types.h"

class Physics
{
public:
    float gravity[3];

    bool readNumber(const std::string &path);
    bool readParameter(const Config &config);
    void initiate();
};

#endif