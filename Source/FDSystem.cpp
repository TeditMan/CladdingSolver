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
    , mu_(inputData_.properties().youngModulus / (2 * (1 + inputData_.properties().poissonRatio)))
    , lambda_(inputData_.properties().youngModulus* inputData_.properties().poissonRatio /
              ((1 + inputData_.properties().poissonRatio) * (1 - 2 * inputData_.properties().poissonRatio)))
{ 
    fillSystem();
}

std::array<std::vector<Triplet>, 2> FDSystem::eq1_left(std::size_t i, std::size_t j)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const double ri{ inputData_.geometry().innerRadius + i * grid_.dr() };

    const double nu{ inputData_.properties().poissonRatio };
    const double mult1{ 2.0 * (1.0 - nu) };
    const double mult2{ 1.0 - 2.0 * nu };

    pushMinorToMajor(discretization_.D_rr_central(i, j), vectorV, mult1);
    pushMinorToMajor(discretization_.D_r_central (i, j), vectorV, mult1 / ri);
    pushMinorToMajor(discretization_.D0          (i, j), vectorV, - mult1 / ri / ri);
    pushMinorToMajor(discretization_.D_zz_central(i, j), vectorV, mult2);
    pushMinorToMajor(discretization_.D_rz_central(i, j), vectorW);

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

    pushMinorToMajor(discretization_.D_rr_central(i, j), vectorW, mult2);
    pushMinorToMajor(discretization_.D_r_central (i, j), vectorW, mult2 / ri);
    pushMinorToMajor(discretization_.D_zz_central(i, j), vectorW, mult1);
    pushMinorToMajor(discretization_.D_z_central (i, j), vectorV, 1 / ri);
    pushMinorToMajor(discretization_.D_rz_central(i, j), vectorV);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::bound1_Ri_left(std::size_t j)
{   
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    pushMinorToMajor(discretization_.D_r_forward(0, j), vectorV, 2 * mu_ + lambda_);
    pushMinorToMajor(discretization_.D0(0, j), vectorV, lambda_ / inputData_.geometry().innerRadius);
    pushMinorToMajor(discretization_.D_z_central(0, j), vectorW, lambda_);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::bound2_Ri_left(std::size_t j)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    pushMinorToMajor(discretization_.D_z_central(0, j), vectorV, mu_);
    pushMinorToMajor(discretization_.D_r_forward(0, j), vectorW, mu_);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::bound1_Ro_left(std::size_t j)
{   
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const std::size_t i{ grid_.Nr() - 1 };

    pushMinorToMajor(discretization_.D_r_backward(i, j), vectorV, 2 * mu_ + lambda_);
    pushMinorToMajor(discretization_.D0(i, j), vectorV, lambda_ / inputData_.geometry().outerRadius);
    pushMinorToMajor(discretization_.D_z_central(i, j), vectorW, lambda_);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::bound2_Ro_left(std::size_t j)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const std::size_t i{ grid_.Nr() - 1 };

    pushMinorToMajor(discretization_.D_z_central(i, j), vectorV, mu_);
    pushMinorToMajor(discretization_.D_r_backward(i, j), vectorW, mu_);

    return { vectorV, vectorW };
}

std::array<std::vector<Triplet>, 2> FDSystem::bound2_H_left(std::size_t i)
{
    std::vector<Triplet> vectorV{};
    std::vector<Triplet> vectorW{};

    const std::size_t j{ grid_.Nz() - 1 };
    const double ri{ inputData_.geometry().innerRadius + i * grid_.dr() };

    pushMinorToMajor(discretization_.D_z_backward(i, j), vectorW, (2 * mu_ + lambda_));
    pushMinorToMajor(discretization_.D_r_central(i, j), vectorV, lambda_);
    pushMinorToMajor(discretization_.D0(i, j), vectorV, lambda_ / ri);

    return { vectorV, vectorW };
}

void FDSystem::pushMinorToMajor(const std::vector<Triplet>& minor, std::vector<Triplet>& major, double mult)
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

void FDSystem::addVectorToA_(const std::vector<Triplet>& Vector, std::size_t row, int type)
{
    for (Triplet tr : Vector)
    {
        const std::size_t column{ 2 * grid_.convertIndices(tr.i, tr.j) + type };
        A_.push_back({ row, column, tr.value });
    }
}

void FDSystem::implementEquations()
{
    for (std::size_t j{ 1 }; j < grid_.Nz() - 1; ++j)
    {
        for (std::size_t i{ 1 }; i < grid_.Nr() - 1; ++i)
        {
            const std::size_t row1{ 2 * grid_.convertIndices(i, j) };
            const std::size_t row2{ row1 + 1 };

            auto [V1, W1] = eq1_left(i, j);
            auto [V2, W2] = eq2_left(i, j);

            addVectorToA_(V1, row1, 0);
            addVectorToA_(W1, row1, 1);
            addVectorToA_(V2, row2, 0);
            addVectorToA_(W2, row2, 1);
        }
    }
}

void FDSystem::implementCorners()
{
    const std::size_t left_bottom_1{ 2 * grid_.convertIndices(0, 0) };
    const std::size_t left_bottom_2{ left_bottom_1 + 1 };
    const std::size_t right_bottom_1{ 2 * grid_.convertIndices(grid_.Nr() - 1, 0) };
    const std::size_t right_bottom_2{ right_bottom_1 + 1 };

    A_.push_back({ right_bottom_1, right_bottom_1, 1 });
    b_[right_bottom_1] = 0.0;
    A_.push_back({ right_bottom_2, right_bottom_2, 1 });
    b_[right_bottom_2] = 0.0;

    A_.push_back({ left_bottom_1, left_bottom_1, 1 });
    b_[left_bottom_1] = 0.0;
    A_.push_back({ left_bottom_2, left_bottom_2, 1 });
    b_[left_bottom_2] = 0.0;

    const std::size_t left_top_1{ 2 * grid_.convertIndices(0, grid_.Nz() - 1) };
    const std::size_t left_top_2{ left_top_1 + 1 };
    const std::size_t right_top_1{ 2 * grid_.convertIndices(grid_.Nr() - 1, grid_.Nz() - 1) };
    const std::size_t right_top_2{ right_top_1 + 1 };

    A_.push_back({ left_top_1, left_top_1, 1 });
    b_[left_top_1] = 0.0;
    A_.push_back({ right_top_1, right_top_1, 1 });
    b_[right_top_1] = 0.0;

    std::vector<Triplet> vectorV_left{};
    std::vector<Triplet> vectorW_left{};
    std::vector<Triplet> vectorV_right{};
    std::vector<Triplet> vectorW_right{};

    const std::size_t j{ grid_.Nz() - 1 };

    pushMinorToMajor(discretization_.D_z_backward(0, j), vectorW_left, (2 * mu_ + lambda_));
    pushMinorToMajor(discretization_.D_r_forward(0, j), vectorV_left, lambda_);
    pushMinorToMajor(discretization_.D0(0, j), vectorV_left, lambda_ / inputData_.geometry().innerRadius);

    pushMinorToMajor(discretization_.D_z_backward(grid_.Nr() - 1, j), vectorW_right, (2 * mu_ + lambda_));
    pushMinorToMajor(discretization_.D_r_backward(grid_.Nr() - 1, j), vectorV_right, lambda_);
    pushMinorToMajor(discretization_.D0(grid_.Nr() - 1, j), vectorV_right, lambda_ / inputData_.geometry().outerRadius);

    addVectorToA_(vectorV_left, left_top_2, 0);
    addVectorToA_(vectorW_left, left_top_2, 1);
    b_[left_top_2] = inputData_.conditions().sigma_zz_H;
    addVectorToA_(vectorV_right, right_top_2, 0);
    addVectorToA_(vectorW_right, right_top_2, 1);
    b_[right_top_2] = inputData_.conditions().sigma_zz_H;
}

void FDSystem::implementBoundaries()
{
    for (std::size_t j{ 1 }; j < grid_.Nz() - 1; ++j)
    {
        const std::size_t row1_Ri{ 2 * grid_.convertIndices(0, j) };
        const std::size_t row2_Ri{ row1_Ri + 1 };

        const std::size_t row1_Ro{ 2 * grid_.convertIndices(grid_.Nr() - 1, j)};
        const std::size_t row2_Ro{ row1_Ro + 1 };


        auto [V1_Ri, W1_Ri] = bound1_Ri_left(j);
        auto [V2_Ri, W2_Ri] = bound2_Ri_left(j);
        auto [V1_Ro, W1_Ro] = bound1_Ro_left(j);
        auto [V2_Ro, W2_Ro] = bound2_Ro_left(j);

        addVectorToA_(V1_Ri, row1_Ri, 0);
        addVectorToA_(W1_Ri, row1_Ri, 1);
        b_[row1_Ri] = inputData_.conditions().sigma_rr_Ri;
        addVectorToA_(V2_Ri, row2_Ri, 0);
        addVectorToA_(W2_Ri, row2_Ri, 1);
        b_[row2_Ri] = inputData_.conditions().sigma_rz_Ri;
        addVectorToA_(V1_Ro, row1_Ro, 0);
        addVectorToA_(W1_Ro, row1_Ro, 1);
        b_[row1_Ro] = inputData_.conditions().sigma_rr_Ro;
        addVectorToA_(V2_Ro, row2_Ro, 0);
        addVectorToA_(W2_Ro, row2_Ro, 1);
        b_[row2_Ro] = inputData_.conditions().sigma_rz_Ro;
    }

    for (std::size_t i{ 1 }; i < grid_.Nr() - 1; ++i)
    {
        const std::size_t row1_0{ 2 * grid_.convertIndices(i, 0) };
        const std::size_t row2_0{ row1_0 + 1 };

        const std::size_t row1_H{ 2 * grid_.convertIndices(i, grid_.Nz() - 1) };
        const std::size_t row2_H{ row1_H + 1 };


        A_.push_back({ row1_0, row1_0, 1 });
        b_[row1_0] = 0.0;
        A_.push_back({ row2_0, row2_0, 1 });
        b_[row2_0] = 0.0;

        A_.push_back({ row1_H, row1_H, 1 });
        b_[row1_H] = 0.0;

        auto [V, W] = bound2_H_left(i);
        addVectorToA_(V, row2_H, 0);
        addVectorToA_(W, row2_H, 1);
        b_[row2_H] = inputData_.conditions().sigma_zz_H;
    }

    implementCorners();
}

void FDSystem::fillSystem()
{
    A_.clear();
    b_.assign(2 * grid_.Nr() * grid_.Nz(), 0.0);

    implementEquations();
    implementBoundaries();
}
