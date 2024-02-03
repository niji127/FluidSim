#include "physics.h"

void Physics::initiate() {}
bool Physics::readNumber(const std::string &path)
{
    return false;
}

bool Physics::readParameter(const Config &config)
{
    {
        float gravity_y = config.Read("gravity", 9.81f);
        gravity[0] = 0.0f;
        gravity[1] = -gravity_y;
        gravity[2] = 0.0f;

        return false;
    }
}