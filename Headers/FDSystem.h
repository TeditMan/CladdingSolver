#pragma once

#include <cstddef>
#include <vector>
#include <array>

#include "Discretization.h"
#include "Grid.h"
#include "InputData.h"
#include "Triplet.h"

class FDSystem
{
public:
    FDSystem(const Discretization&, const Grid&, const InputData&);

    void implementEquations();
    void implementBoundaries();
    void fillSystem();

private:
    const Discretization& discretization_;
    const Grid& grid_;
    const InputData& inputData_;
    
    const double mu_;
    const double lambda_;
    
    std::vector<Triplet> A_{};
    std::vector<double>  b_{};

    std::array<std::vector<Triplet>, 2> eq1_left(std::size_t, std::size_t);
    std::array<std::vector<Triplet>, 2> eq2_left(std::size_t, std::size_t);

    std::array<std::vector<Triplet>, 2> bound1_Ri_left(std::size_t);
    std::array<std::vector<Triplet>, 2> bound2_Ri_left(std::size_t);

    std::array<std::vector<Triplet>, 2> bound1_Ro_left(std::size_t);
    std::array<std::vector<Triplet>, 2> bound2_Ro_left(std::size_t);

    void addVectorToA_(const std::vector<Triplet>&, std::size_t, int);
    void pushMinorToMajor(const std::vector<Triplet>&, std::vector<Triplet>&, double mult = 1);

};
