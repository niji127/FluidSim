#include "output.cuh"
#include "thrust/host_vector.h"
#include <windows.h>
#include <boost/iostreams/device/mapped_file.hpp>

void Output::writeFluidVelocity(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number;
    thrust::host_vector<glm::vec3> fluid_position = sim->fluid_list.position;
    thrust::host_vector<glm::vec3> fluid_velocity = sim->fluid_list.velocity;

    std::string file_path = path + "\\fluid_velocity-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "Velocity vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(fluid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(fluid_num) + "\n";
    str2 += "VECTORS velocity float\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_position[i].x) + " " + std::to_string(fluid_position[i].y) + " " + std::to_string(fluid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_velocity[i].x) + " " + std::to_string(fluid_velocity[i].y) + " " + std::to_string(fluid_velocity[i].z) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeFluidDensity(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number;
    thrust::host_vector<glm::vec3> fluid_position = sim->fluid_list.position;
    thrust::host_vector<float> fluid_density = sim->fluid_list.density;

    std::string file_path = path + "\\fluid_density-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "scalar vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(fluid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(fluid_num) + "\n";
    str2 += "SCALARS density float 1\n";
    str2 += "LOOKUP_TABLE default\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_position[i].x) + " " + std::to_string(fluid_position[i].y) + " " + std::to_string(fluid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_density[i]) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeFluidPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int fluid_num = Sim::parameter.fluid.number;
    thrust::host_vector<glm::vec3> fluid_position = sim->fluid_list.position;
    thrust::host_vector<float> fluid_pressure = sim->fluid_list.color;

    std::string file_path = path + "\\fluid_pressure-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "scalar vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(fluid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(fluid_num) + "\n";
    str2 += "SCALARS pressure float 1\n";
    str2 += "LOOKUP_TABLE default\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_position[i].x) + " " + std::to_string(fluid_position[i].y) + " " + std::to_string(fluid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < fluid_num; ++i)
        str += std::to_string(fluid_pressure[i]) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeSolidVelocity(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<glm::vec3> solid_position = sim->solid_list.position;
    thrust::host_vector<glm::vec3> solid_velocity = sim->solid_list.velocity;

    std::string file_path = path + "\\solid_velocity-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "Velocity vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(solid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(solid_num) + "\n";
    str2 += "VECTORS velocity float\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_position[i].x) + " " + std::to_string(solid_position[i].y) + " " + std::to_string(solid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_velocity[i].x) + " " + std::to_string(solid_velocity[i].y) + " " + std::to_string(solid_velocity[i].z) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeSolidPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<glm::vec3> solid_position = sim->solid_list.position;
    thrust::host_vector<float> solid_pressure = sim->solid_list.pressure;

    std::string file_path = path + "\\solid_pressure-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "scalar vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(solid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(solid_num) + "\n";
    str2 += "SCALARS pressure float 1\n";
    str2 += "LOOKUP_TABLE default\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_position[i].x) + " " + std::to_string(solid_position[i].y) + " " + std::to_string(solid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_pressure[i]) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeSolidDensity(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    float density = Sim::parameter.solid.reference_density;
    thrust::host_vector<glm::vec3> solid_position = sim->solid_list.position;

    std::string file_path = path + "\\solid_density-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "scalar vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(solid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(solid_num) + "\n";
    str2 += "SCALARS density float 1\n";
    str2 += "LOOKUP_TABLE default\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_position[i].x) + " " + std::to_string(solid_position[i].y) + " " + std::to_string(solid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(density) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeSolidNormal(Sim *sim, const std::string &path, SimTime *time)
{
    int solid_num = Sim::parameter.solid.number;
    thrust::host_vector<glm::vec3> solid_position = sim->solid_list.position;
    thrust::host_vector<glm::vec3> solid_normal = sim->solid_list.normal;

    std::string file_path = path + "\\solid_normal-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "Velocity vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(solid_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(solid_num) + "\n";
    str2 += "VECTORS normal float\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_position[i].x) + " " + std::to_string(solid_position[i].y) + " " + std::to_string(solid_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < solid_num; ++i)
        str += std::to_string(solid_normal[i].x) + " " + std::to_string(solid_normal[i].y) + " " + std::to_string(solid_normal[i].z) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeVirtNormal(Sim *sim, const std::string &path, SimTime *time)
{
    int virt_num = Sim::parameter.virt.number;
    thrust::host_vector<glm::vec3> virt_position = sim->virt_list.position;
    thrust::host_vector<glm::vec3> virt_normal = sim->virt_list.normal;

    std::string file_path = path + "\\virt_normal-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "Velocity vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(virt_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(virt_num) + "\n";
    str2 += "VECTORS normal float\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < virt_num; ++i)
        str += std::to_string(virt_position[i].x) + " " + std::to_string(virt_position[i].y) + " " + std::to_string(virt_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < virt_num; ++i)
        str += std::to_string(virt_normal[i].x) + " " + std::to_string(virt_normal[i].y) + " " + std::to_string(virt_normal[i].z) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::writeVirtPressure(Sim *sim, const std::string &path, SimTime *time)
{
    int virt_num = Sim::parameter.virt.number;
    thrust::host_vector<glm::vec3> virt_position = sim->virt_list.position;
    thrust::host_vector<float> virt_pressure = sim->virt_list.pressure;

    std::string file_path = path + "\\virt_pressure-" + std::to_string(time->file_step) + ".vtk";
    std::size_t file_size = 500 * 1024 * 1024;

    std::string str1 = "# vtk DataFile Version 2.0\n";
    str1 += "scalar vector file\n";
    str1 += "ASCII\n";
    str1 += "DATASET POLYDATA\n";
    str1 += "POINTS " + std::to_string(virt_num) + " float\n";

    std::string str2 = "POINT_DATA " + std::to_string(virt_num) + "\n";
    str2 += "SCALARS pressure float 1\n";
    str2 += "LOOKUP_TABLE default\n";

    std::string str = str1;
    str.reserve(file_size);
    for (int i = 0; i < virt_num; ++i)
        str += std::to_string(virt_position[i].x) + " " + std::to_string(virt_position[i].y) + " " + std::to_string(virt_position[i].z) + "\n";
    str += str2;
    for (int i = 0; i < virt_num; ++i)
        str += std::to_string(virt_pressure[i]) + "\n";

    std::fstream output;
    output.open(file_path, std::ios::out);
    output << str;
    output.close();
}

void Output::outputData(Sim *sim)
{
    SimTime *time = Sim::parameter.getTime();
    if (!(time->isOutputData()))
        return;
    clock_t start = clock();
    std::cout << "Outputing Data..." << std::endl;
    std::cout << "--time=" << time->current_time - time->dt << "s step=" << time->i - 1 << "--" << std::endl;
    std::string file_path = "..\\..\\out";

    if (Sim::parameter.hasFluid())
    {
        std::string fluid_file_path = file_path + "\\fluid";
        CreateDirectory(fluid_file_path.c_str(), NULL);
        writeFluidVelocity(sim, fluid_file_path, time);
        writeFluidDensity(sim, fluid_file_path, time);
        writeFluidPressure(sim, fluid_file_path, time);
    }
    if (Sim::parameter.hasSolid())
    {
        std::string solid_file_path = file_path + "\\solid";
        CreateDirectory(solid_file_path.c_str(), NULL);
        writeSolidVelocity(sim, solid_file_path, time);
        writeSolidPressure(sim, solid_file_path, time);
        writeSolidDensity(sim, solid_file_path, time);
        writeSolidNormal(sim, solid_file_path, time);
    }
    if (Sim::parameter.hasVirtual())
    {
        std::string virt_file_path = file_path + "\\virtual";
        CreateDirectory(virt_file_path.c_str(), NULL);
        writeVirtNormal(sim, virt_file_path, time);
        writeVirtPressure(sim, virt_file_path, time);
    }
    clock_t end = clock();
    std::cout << "write time=" << end - start << "ms" << std::endl;
}