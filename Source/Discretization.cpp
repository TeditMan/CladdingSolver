#include <vector>

#include "Discretization.h"
#include "Grid.h"
#include "Triplet.h"

std::vector<Triplet> CentralD::D0(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j, 1.0 }
    };
}

std::vector<Triplet> CentralD::D_r(std::size_t i, std::size_t j) const
{
    return 
    { 
        {i + 1, j, (+1.0 / (2.0 * grid_.dr())) },
        {i - 1, j, (-1.0 / (2.0 * grid_.dr())) } 
    };
}

std::vector<Triplet> CentralD::D_z(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 1, (+1.0 / (2.0 * grid_.dz())) },
        { i, j - 1, (-1.0 / (2.0 * grid_.dz())) }
    };
}

std::vector<Triplet> CentralD::D_rr(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 1, j, (+1.0 / (grid_.dr() * grid_.dr())) },
        { i    , j, (-2.0 / (grid_.dr() * grid_.dr())) },
        { i - 1, j, (+1.0 / (grid_.dr() * grid_.dr())) }
    };
}

std::vector<Triplet> CentralD::D_zz(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 1, (+1.0 / (grid_.dz() * grid_.dz())) },
        { i, j    , (-2.0 / (grid_.dz() * grid_.dz())) },
        { i, j - 1, (+1.0 / (grid_.dz() * grid_.dz())) }
    };
}

std::vector<Triplet> CentralD::D_rz(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 1, j + 1, (+1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i + 1, j - 1, (-1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i - 1, j + 1, (-1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i - 1, j - 1, (+1.0 / (4.0 * grid_.dr() * grid_.dz())) }
    };
}
