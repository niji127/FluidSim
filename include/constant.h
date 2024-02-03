#include "parameter.h"
#include "particle.h"

struct PhysicsConstant
{
    float gravity[3];
};

struct FluidConstant
{
    int number;
    float sound_speed, reference_density;
    float reference_density_inverse;
    float gamma, viscosity; // dynamic viscosity
    float gamma_inv;
    float coefficient_p2rho;
    float coefficient_rho2p;
    float min_pressure;
    bool is_fluid_viscid;
};

struct SolidConstant
{
    int number;
    float density, modulus, poisson_ratio;
    float artificial_viscocity[2];
};

struct VirtConstant
{
    int number;
};

struct Domain
{
    float pos_min[3], pos_max[3];
    float boundary_min[3], boundary_max[3];
    int block_number[3];
    float interval[3];
    float interval_inv[3];
    int block_number_total;
};

struct TimeConstant
{
    float dt, current_time;
    int istart, iend, i;
    int file_step;
    int solid_sub_step, shift_sub_step;
    int result_per_step;
}

class Kernel
{
public:
    float smoothing_length;
    float particle_diameter;
    float particle_volume;
    float diameter_square;

    int block_size;
    int pair_size;
    int thread_num;

    float smoothing_length_inverse;

    float smoothing_length_inverse_square;

    float impact_distance_square;

    float kernel_coefficient_1;
    float kernel_coefficient_2;

    float kernel_differential;

    bool readParameter(const Config &config);
    void initiate();
    bool readNumber(const std::string &path);
};

class CudaConstant
{
public:
    float smoothing_length;
    float particle_diameter;
    float particle_volume;
    float diameter_square;

    int block_size;
    int fluid_pair_size;
    int solid_pair_size;
    int thread_num;

    float impact_distance_square;

    float kernel_coefficient_1;
    float kernel_coefficient_2;

    float kernel_differential;
    float coefficient_3to2;

    Fluid_Particle *fluid_particle;
    Solid_Particle *solid_particle;
    Virtual_Particle *virtual_particle;

    float pos_min[3], pos_max[3];
    float gravity[3];

    void copyFromPhysics(Physics *physics)
    {
        for (int i = 0; i < 3; i++)
            this->gravity[i] = physics->gravity[i];
    }

    void copyFromDomain(Domain *domain)
    {
        for (int i = 0; i < 3; i++)
        {
            this->pos_min[i] = domain->boundary_min[i];
            this->pos_max[i] = domain->boundary_max[i];
        }
    }

    void copyFromKernel(Kernel *kernel)
    {
        this->smoothing_length = kernel->smoothing_length;
        this->particle_diameter= kernel->particle_diameter;
        this->particle_volume  = kernel->particle_volume;
        this->diameter_square  = kernel->diameter_square;

        this->block_size = kernel->block_size;
        this->pair_size = kernel->pair_size;
        this->thread_num = kernel->thread_num;

        this->smoothing_length_inverse = kernel->smoothing_length_inverse;

        this->smoothing_length_inverse_square = kernel->smoothing_length_inverse_square;

        this->impact_distance_square = kernel->impact_distance_square;

        this->kernel_coefficient_1 = kernel->kernel_coefficient_1;
        this->kernel_coefficient_2 = kernel->kernel_coefficient_2;

        this->kernel_differential = kernel->kernel_differential;
        this->coefficient_3to2 = kernel->coefficient_3to2;
    }
};