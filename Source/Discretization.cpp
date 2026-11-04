#include <vector>

#include "Discretization.h"
#include "Grid.h"
#include "Triplet.h"

std::vector<Triplet> Discretization::D0(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j, 1.0 }
    };
}



std::vector<Triplet> Discretization::D_r_central(std::size_t i, std::size_t j) const
{
    return 
    { 
        {i + 1, j, (+1.0 / (2.0 * grid_.dr())) },
        {i - 1, j, (-1.0 / (2.0 * grid_.dr())) } 
    };
}

std::vector<Triplet> Discretization::D_z_central(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 1, (+1.0 / (2.0 * grid_.dz())) },
        { i, j - 1, (-1.0 / (2.0 * grid_.dz())) }
    };
}

std::vector<Triplet> Discretization::D_rr_central(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 1, j, (+1.0 / (grid_.dr() * grid_.dr())) },
        { i    , j, (-2.0 / (grid_.dr() * grid_.dr())) },
        { i - 1, j, (+1.0 / (grid_.dr() * grid_.dr())) }
    };
}

std::vector<Triplet> Discretization::D_zz_central(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 1, (+1.0 / (grid_.dz() * grid_.dz())) },
        { i, j    , (-2.0 / (grid_.dz() * grid_.dz())) },
        { i, j - 1, (+1.0 / (grid_.dz() * grid_.dz())) }
    };
}

std::vector<Triplet> Discretization::D_rz_central(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 1, j + 1, (+1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i + 1, j - 1, (-1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i - 1, j + 1, (-1.0 / (4.0 * grid_.dr() * grid_.dz())) },
        { i - 1, j - 1, (+1.0 / (4.0 * grid_.dr() * grid_.dz())) }
    };
}



std::vector<Triplet> Discretization::D_r_forward(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 2, j, (-1.0 / grid_.dr()) },
        { i + 1, j, (+4.0 / grid_.dr()) },
        { i    , j, (-3.0 / grid_.dr()) }
    };
}

std::vector<Triplet> Discretization::D_z_forward(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 2, (-1.0 / grid_.dz()) },
        { i, j + 1, (+4.0 / grid_.dz()) },
        { i, j    , (-3.0 / grid_.dz()) }
    };
}

std::vector<Triplet> Discretization::D_r_backward(std::size_t i, std::size_t j) const
{
    return
    {
        { i + 2, j, (+1.0 / grid_.dr()) },
        { i + 1, j, (-4.0 / grid_.dr()) },
        { i    , j, (+3.0 / grid_.dr()) }
    };
}

std::vector<Triplet> Discretization::D_z_backward(std::size_t i, std::size_t j) const
{
    return
    {
        { i, j + 2, (+1.0 / grid_.dz()) },
        { i, j + 1, (-4.0 / grid_.dz()) },
        { i, j    , (+3.0 / grid_.dz()) }
    };
}
