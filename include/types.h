#ifndef _TYPES_H
#define _TYPES_H

enum ParticleType
{
    FLUID,
    SOLID,
    VIRT
};

enum SolidType
{
    UNCONSTRAINED_SOLID,
    FIXED_SOLID
};

enum VirtType
{
    FIXED_VIRT,
    MOVING_VIRT
};

#endif