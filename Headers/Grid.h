#pragma once

#include <cstddef>

#include "InputData.h"

struct NodesData
{
    std::size_t Nodes_r_{};
    std::size_t Nodes_z_{};
};

class Grid
{
public:
    Grid(const InputData&, const NodesData&);

    std::size_t convertIndices(std::size_t, std::size_t) const;

    double dr() const { return dr_; }
    double dz() const { return dz_; }
    std::size_t Nr() const { return Nr_; }
    std::size_t Nz() const { return Nz_; }

private:
    std::size_t Nr_{};
    std::size_t Nz_{};

    double dr_{};
    double dz_{};

    double Ri_{};
    double Ro_{};
    double H_{};
};
