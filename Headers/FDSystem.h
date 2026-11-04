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

    void fillMatrix();

private:
    std::vector<Triplet> A_{};
    std::vector<double> b_{};

    const Discretization& discretization_;
    const Grid& grid_;
    const InputData& inputData_;

    std::array<std::vector<Triplet>, 2> eq1_left(std::size_t, std::size_t);
    std::array<std::vector<Triplet>, 2> eq2_left(std::size_t, std::size_t);

    void pushMinorToMajor(const std::vector<Triplet>&, std::vector<Triplet>&, double mult = 1);

};
