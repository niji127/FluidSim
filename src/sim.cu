#include "sim.cuh"
#include "convert.h"

#include "host_data.cuh"
#include "thrust/host_vector.h"

#include <string>

Parameter Sim::parameter;
Data_Pointer Sim::data_pointer;

void Sim::showMemoryUsed()
{
    // particles
    long long fluid_size = fluid_list.getMemoryUsed();
    long long solid_size = solid_list.getMemoryUsed();
    long long virt_size = virt_list.getMemoryUsed();

    // cells
    thrust::host_vector<Cell> host_cell;
    host_cell.push_back(cell_list[0]);
    long long cell_size = host_cell[0].getMemoryUsed();
    cell_size *= cell_list.size();

    long long parameter_size = sizeof(Parameter);
    long long sim_size = sizeof(Sim);
    long long miscellaneous_size = parameter_size + sim_size;

    long long all_size = fluid_size + solid_size + virt_size + cell_size + miscellaneous_size;

    long long digit = 2;
    std::cout << "\nAll memory: " << convertUnit(all_size, digit) << std::endl;
    std::cout << "fluid particle: " << convertUnit(fluid_size, digit) << std::endl;
    std::cout << "solid particle: " << convertUnit(solid_size, digit) << std::endl;
    std::cout << "virtual particle: " << convertUnit(virt_size, digit) << std::endl;
    std::cout << "cell: " << convertUnit(cell_size, digit) << std::endl;
    std::cout << "other: " << convertUnit(miscellaneous_size, digit) << std::endl;
}

void Sim::readParticleData()
{
    std::cout << "\nReading Particles..." << std::endl;
    bool error_flag(false);
    std::string path = "..\\..\\data";

    clock_t start, end;
    start = clock();

    error_flag = readFluidData(path);
    error_flag = readSolidData(path) || error_flag;
    error_flag = readVirtData(path) || error_flag;
    end = clock();
    printf("time=%dms\n", end - start);

    if (error_flag)
        throw std::invalid_argument("FAILED TO READ PARTICLES DATA");
}

bool Sim::readFluidData(const std::string &path)
{
    int fluid_num = parameter.fluid.number;
    if (fluid_num == 0)
    {
        std::cout << "\nno fluid particle" << std::endl;
        return false;
    }

    if (fluid_num < 0)
    {
        std::cout << "fluid particle number should >= 0" << std::endl;
        return true;
    }

    std::ifstream infile;
    infile.open(path + "\\fluid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read fluid particles data" << std::endl;
        return true;
    }

    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        fluid_host_list.addList(str, parameter);
    infile.close();

    return false;
}

bool Sim::readSolidData(const std::string &path)
{
    int solid_num = parameter.solid.number;
    if (solid_num == 0)
    {
        std::cout << "\nno solid particle" << std::endl;
        return false;
    }
    if (solid_num < 0)
    {
        std::cout << "solid particle number should >= 0" << std::endl;
        return true;
    }
    std::ifstream infile;
    infile.open(path + "\\solid_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read solid particles data" << std::endl;
        return true;
    }

    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        solid_host_list.addList(str, parameter);
    infile.close();

    return false;
}

bool Sim::readVirtData(const std::string &path)
{
    int virt_num = parameter.virt.number;
    if (virt_num == 0)
    {
        std::cout << "\nno virtual particle" << std::endl;
        return false;
    }

    if (virt_num < 0)
    {
        std::cout << "virtual particle number should >= 0" << std::endl;
        return true;
    }

    std::ifstream infile;
    infile.open(path + "\\virt_particle.dat", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "failed to read virtual particles data" << std::endl;
        return true;
    }

    std::string str;
    std::getline(infile, str);
    while (std::getline(infile, str))
        virt_host_list.addList(str, parameter);
    infile.close();

    return false;
}

void Sim::createCell()
{
    int cell_number = parameter.domain.cell_number_total;
    int(&cell_n)[3] = parameter.domain.cell_number;
    cell_host_list.reserve(cell_number);

    for (int i = 0; i < cell_n[0]; i++)
    {
        for (int j = 0; j < cell_n[1]; j++)
        {
            for (int k = 0; k < cell_n[2]; k++)
            {
                Cell temp_cell;
                int id = i * cell_n[1] * cell_n[2] + j * cell_n[2] + k;
                if (i != 0)
                {
                    if ((j != 0) && (k != 0))
                        temp_cell.neighbor_cell_id[0] = id - cell_n[1] * cell_n[2] - cell_n[2] - 1;
                    if ((j != 0) && (k != cell_n[2] - 1))
                        temp_cell.neighbor_cell_id[1] = id - cell_n[1] * cell_n[2] - cell_n[2] + 1;
                    if ((j != cell_n[1] - 1) && (k != 0))
                        temp_cell.neighbor_cell_id[2] = id - cell_n[1] * cell_n[2] + cell_n[2] - 1;
                    if ((j != cell_n[1] - 1) && (k != cell_n[2] - 1))
                        temp_cell.neighbor_cell_id[3] = id - cell_n[1] * cell_n[2] + cell_n[2] + 1;

                    if (k != 0)
                        temp_cell.neighbor_cell_id[4] = id - cell_n[1] * cell_n[2] - 1;
                    if (k != cell_n[2] - 1)
                        temp_cell.neighbor_cell_id[5] = id - cell_n[1] * cell_n[2] + 1;
                    if (j != 0)
                        temp_cell.neighbor_cell_id[6] = id - cell_n[1] * cell_n[2] - cell_n[2];
                    if (j != cell_n[1] - 1)
                        temp_cell.neighbor_cell_id[7] = id - cell_n[1] * cell_n[2] + cell_n[2];

                    temp_cell.neighbor_cell_id[8] = id - cell_n[1] * cell_n[2];
                }

                if ((j != 0) && (k != 0))
                    temp_cell.neighbor_cell_id[9] = id - cell_n[2] - 1;
                if ((j != 0) && (k != cell_n[2] - 1))
                    temp_cell.neighbor_cell_id[10] = id - cell_n[2] + 1;
                if ((j != cell_n[1] - 1) && (k != 0))
                    temp_cell.neighbor_cell_id[11] = id + cell_n[2] - 1;
                if ((j != cell_n[1] - 1) && (k != cell_n[2] - 1))
                    temp_cell.neighbor_cell_id[12] = id + cell_n[2] + 1;

                if (k != 0)
                    temp_cell.neighbor_cell_id[13] = id - 1;
                if (k != cell_n[2] - 1)
                    temp_cell.neighbor_cell_id[14] = id + 1;
                if (j != 0)
                    temp_cell.neighbor_cell_id[15] = id - cell_n[2];
                if (j != cell_n[1] - 1)
                    temp_cell.neighbor_cell_id[16] = id + cell_n[2];

                temp_cell.neighbor_cell_id[17] = id;

                if (i != cell_n[0] - 1)
                {
                    if ((j != 0) && (k != 0))
                        temp_cell.neighbor_cell_id[18] = id + cell_n[1] * cell_n[2] - cell_n[2] - 1;
                    if ((j != 0) && (k != cell_n[2] - 1))
                        temp_cell.neighbor_cell_id[19] = id + cell_n[1] * cell_n[2] - cell_n[2] + 1;
                    if ((j != cell_n[1] - 1) && (k != 0))
                        temp_cell.neighbor_cell_id[20] = id + cell_n[1] * cell_n[2] + cell_n[2] - 1;
                    if ((j != cell_n[1] - 1) && (k != cell_n[2] - 1))
                        temp_cell.neighbor_cell_id[21] = id + cell_n[1] * cell_n[2] + cell_n[2] + 1;

                    if (k != 0)
                        temp_cell.neighbor_cell_id[22] = id + cell_n[1] * cell_n[2] - 1;
                    if (k != cell_n[2] - 1)
                        temp_cell.neighbor_cell_id[23] = id + cell_n[1] * cell_n[2] + 1;
                    if (j != 0)
                        temp_cell.neighbor_cell_id[24] = id + cell_n[1] * cell_n[2] - cell_n[2];
                    if (j != cell_n[1] - 1)
                        temp_cell.neighbor_cell_id[25] = id + cell_n[1] * cell_n[2] + cell_n[2];

                    temp_cell.neighbor_cell_id[26] = id + cell_n[1] * cell_n[2];
                }
                cell_host_list.push_back(temp_cell);
            }
        }
    }
}

int Sim::getCellID(const glm::vec3 &pos)
{
    float xmin = parameter.domain.cell_min[0];
    float ymin = parameter.domain.cell_min[1];
    float zmin = parameter.domain.cell_min[2];
    float inter_inv_x = parameter.domain.interval_inv[0];
    float inter_inv_y = parameter.domain.interval_inv[1];
    float inter_inv_z = parameter.domain.interval_inv[2];
    int block_ny = parameter.domain.cell_number[1];
    int block_nz = parameter.domain.cell_number[2];

    int nx = (int)((pos.x - xmin) * inter_inv_x);
    int ny = (int)((pos.y - ymin) * inter_inv_y);
    int nz = (int)((pos.z - zmin) * inter_inv_z);

    return nx * block_ny * block_nz + ny * block_nz + nz;
}

void Sim::fillCell()
{
    int fluid_num = parameter.fluid.number;
    for (int i = 0; i < fluid_num; i++)
    {
        int cell_id = getCellID(fluid_host_list.position[i]);
        if (cell_host_list[cell_id].fluid_start_id == -1)
            cell_host_list[cell_id].fluid_start_id = i;
        else
        {
            int current_id = cell_host_list[cell_id].fluid_start_id;
            while (fluid_host_list.next_id[current_id] != -1)
                current_id = fluid_host_list.next_id[current_id];
            fluid_host_list.next_id[current_id] = i;
        }
    }

    int solid_num = parameter.solid.number;
    for (int i = 0; i < solid_num; i++)
    {
        int cell_id = getCellID(solid_host_list.position[i]);
        if (cell_host_list[cell_id].solid_start_id == -1)
            cell_host_list[cell_id].solid_start_id = i;
        else
        {
            int current_id = cell_host_list[cell_id].solid_start_id;
            while (solid_host_list.next_id[current_id] != -1)
                current_id = solid_host_list.next_id[current_id];
            solid_host_list.next_id[current_id] = i;
        }
    }

    int virt_num = parameter.virt.number;
    for (int i = 0; i < virt_num; i++)
    {
        int cell_id = getCellID(virt_host_list.position[i]);
        if (cell_host_list[cell_id].virt_start_id == -1)
            cell_host_list[cell_id].virt_start_id = i;
        else
        {
            int current_id = cell_host_list[cell_id].virt_start_id;
            while (virt_host_list.next_id[current_id] != -1)
                current_id = virt_host_list.next_id[current_id];
            virt_host_list.next_id[current_id] = i;
        }
    }
}

void Sim::copyData()
{
    fluid_list.initiate(parameter);
    solid_list.initiate(parameter);
    virt_list.initiate(parameter);

    fluid_list.copyData(fluid_host_list);
    solid_list.copyData(solid_host_list);
    virt_list.copyData(virt_host_list);
    cell_list = cell_host_list;
}
void Sim::setPointer()
{
    data_pointer.fluidInitiate(fluid_list);
    data_pointer.solidInitiate(solid_list);
    data_pointer.virtInitiate(virt_list);
    data_pointer.cellInitiate(cell_list);
}

void Sim::clearHost()
{
    fluid_host_list.clear();
    solid_host_list.clear();
    virt_host_list.clear();
    cell_host_list.clear();
    cell_host_list.shrink_to_fit();
}