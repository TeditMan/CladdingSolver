#pragma once

#include <cstddef>
#include <vector>

#include "Grid.h"
#include "Triplet.h"

class Discretization
{
public:
    explicit Discretization(const Grid& grid)
        : grid_(grid)
    { }

    std::vector<Triplet> D0           (std::size_t, std::size_t) const;

    std::vector<Triplet> D_r_central (std::size_t, std::size_t) const;
    std::vector<Triplet> D_z_central (std::size_t, std::size_t) const;
    std::vector<Triplet> D_rr_central(std::size_t, std::size_t) const;
    std::vector<Triplet> D_zz_central(std::size_t, std::size_t) const;
    std::vector<Triplet> D_rz_central(std::size_t, std::size_t) const;

    std::vector<Triplet> D_r_forward (std::size_t, std::size_t) const;
    std::vector<Triplet> D_z_forward (std::size_t, std::size_t) const;
    std::vector<Triplet> D_r_backward(std::size_t, std::size_t) const;
    std::vector<Triplet> D_z_backward(std::size_t, std::size_t) const;

private:
    const Grid& grid_;

};
