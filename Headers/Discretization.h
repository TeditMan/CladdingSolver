#pragma once

#include <cstddef>
#include <vector>

#include "Grid.h"
#include "Triplet.h"

class Discretization
{
public:
    Discretization(const Grid& grid)
        : grid_(grid)
    { }

    virtual ~Discretization() = default;

    virtual std::vector<Triplet> D0  (std::size_t, std::size_t) const = 0;
    virtual std::vector<Triplet> D_r (std::size_t, std::size_t) const = 0;
    virtual std::vector<Triplet> D_z (std::size_t, std::size_t) const = 0;
    virtual std::vector<Triplet> D_rr(std::size_t, std::size_t) const = 0;
    virtual std::vector<Triplet> D_zz(std::size_t, std::size_t) const = 0;
    virtual std::vector<Triplet> D_rz(std::size_t, std::size_t) const = 0;

protected:
    const Grid& grid_;
};

class CentralD : public Discretization
{ 
public:
    explicit CentralD(const Grid& grid)
        : Discretization(grid)
    { }

    virtual std::vector<Triplet> D0  (std::size_t, std::size_t) const override;
    virtual std::vector<Triplet> D_r (std::size_t, std::size_t) const override;
    virtual std::vector<Triplet> D_z (std::size_t, std::size_t) const override;
    virtual std::vector<Triplet> D_rr(std::size_t, std::size_t) const override;
    virtual std::vector<Triplet> D_zz(std::size_t, std::size_t) const override;
    virtual std::vector<Triplet> D_rz(std::size_t, std::size_t) const override;
};
