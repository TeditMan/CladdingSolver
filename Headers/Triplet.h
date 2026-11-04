#pragma once

#include <stdexcept>
#include <cstddef>

struct Triplet
{
    std::size_t i{};
    std::size_t j{};
    double value{};
};

inline Triplet operator+(const Triplet& a, const Triplet& b)
{
    if (a.i != b.i || a.j != b.j)
    {
        throw std::invalid_argument("Cannot add triplets with different indices");
    }

    return { a.i, a.j, a.value + b.value };
}

inline Triplet operator-(const Triplet& a, const Triplet& b)
{
    if (a.i != b.i || a.j != b.j)
    {
        throw std::invalid_argument("Cannot subtract triplets with different indices");
    }

    return { a.i, a.j, a.value - b.value };
}

inline Triplet operator*(double number, const Triplet& t)
{
    return { t.i, t.j, number * t.value };
}

inline Triplet operator*(const Triplet& t, double number)
{
    return number * t;
}

inline Triplet operator/(const Triplet& t, double number)
{
    if (number == 0.0)
    {
        throw std::invalid_argument("Division by zero");
    }

    return { t.i, t.j, t.value / number };
}
