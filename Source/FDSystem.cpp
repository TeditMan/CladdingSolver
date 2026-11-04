#include <cstddef>
#include <vector>
#include <array>


#include "FDSystem.h"
#include "Grid.h"
#include "Discretization.h"
#include "Triplet.h"
#include "InputData.h"


FDSystem::FDSystem(const Discretization& d, const Grid& g, const InputData& i)
    : discretization_(d)
    , grid_(g)
    , inputData_(i)
{ }

std::array<std::vector<Triplet>, 2> FDSystem::eq1_left(std::size_t i, std::size_t j)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const double ri{ inputData_.geometry().innerRadius + i * grid_.dr() };

    const double nu{ inputData_.properties().poissonRatio };
    const double mult1{ 2.0 * (1.0 - nu) };
    const double mult2{ 1.0 - 2.0 * nu };

    pushMinorToMajor(discretization_.D_rr(i, j), vectorV, mult1);
    pushMinorToMajor(discretization_.D_r (i, j), vectorV, mult1 / ri);
    pushMinorToMajor(discretization_.D0  (i, j), vectorV, - mult1 / ri / ri);
    pushMinorToMajor(discretization_.D_zz(i, j), vectorV, mult2);
    pushMinorToMajor(discretization_.D_rz(i, j), vectorW);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::eq2_left(std::size_t i, std::size_t j)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const double ri{ inputData_.geometry().innerRadius + i * grid_.dr() };

    const double nu{ inputData_.properties().poissonRatio };
    const double mult1{ 2.0 * (1.0 - nu) };
    const double mult2{ 1.0 - 2.0 * nu };

    pushMinorToMajor(discretization_.D_rr(i, j), vectorW, mult2);
    pushMinorToMajor(discretization_.D_r (i, j), vectorW, mult2 / ri);
    pushMinorToMajor(discretization_.D_zz(i, j), vectorW, mult1);
    pushMinorToMajor(discretization_.D_z (i, j), vectorV, 1 / ri);
    pushMinorToMajor(discretization_.D_rz(i, j), vectorV);

    return { vectorV, vectorW };
}

void FDSystem::pushMinorToMajor(std::vector<Triplet> minor, std::vector<Triplet>& major, double mult)
{
    for (Triplet cand : minor)
    {
        cand = cand * mult;
        bool found{ false };
        for (auto& old : major)
        {
            if (cand.i == old.i && cand.j == old.j)
            {
                old = old + cand;
                found = true;
                break;
            }
        }

        if (!found)
        {
            major.push_back(cand);
        }
    }
}

void FDSystem::fillMatrix()
{
    A_.clear();
    
    for (std::size_t j{ 1 }; j < grid_.Nz() - 1; ++j)
    {
        for (std::size_t i{ 1 }; i < grid_.Nr() - 1; ++i)
        {
            const std::size_t row1{ 2 * grid_.convertIndices(i, j) };
            const std::size_t row2{ row1 + 1 };

            auto [V1, W1] = eq1_left(i, j);
            auto [V2, W2] = eq2_left(i, j);

            for (Triplet tr : V1)
            {
                const std::size_t column{ 2 * grid_.convertIndices(tr.i, tr.j) };
                A_.push_back({ row1, column, tr.value });
            }

            for (Triplet tr : W1)
            {
                const std::size_t column{ 2 * grid_.convertIndices(tr.i, tr.j) + 1 };
                A_.push_back({ row1, column, tr.value });
            }

            for (Triplet tr : V2)
            {
                const std::size_t column{ 2 * grid_.convertIndices(tr.i, tr.j) };
                A_.push_back({ row2, column, tr.value });
            }

            for (Triplet tr : W2)
            {
                const std::size_t column{ 2 * grid_.convertIndices(tr.i, tr.j) + 1 };
                A_.push_back({ row2, column, tr.value });
            }
        }
    }
}
