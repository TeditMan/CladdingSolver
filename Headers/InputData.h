#pragma once

struct Pressure
{
    double Pgap{};
    double Pcool{};
};

struct Geometry
{
    double innerRadius{};
    double outerRadius{};
    double height{};
};

struct MaterialProperties
{
    double youngModulus{};
    double poissonRatio{};
};

struct BoundaryConditions
{
    static constexpr double v_z0{ 0.0 };
    static constexpr double w_z0{ 0.0 };
    static constexpr double v_zH{ 0.0 };

    static constexpr double sigma_rz_Ri{ 0.0 };
    static constexpr double sigma_rz_Ro{ 0.0 };

    double sigma_rr_Ri{};
    double sigma_rr_Ro{};
    double sigma_zz_H{};
};

class InputData
{
public:
    InputData(Pressure, Geometry, MaterialProperties);

    Geometry           geometry()   const { return geometry_; }
    MaterialProperties properties() const { return properties_; }
    BoundaryConditions conditions() const { return conditions_; }

private:
    Pressure pressure_{};
    Geometry geometry_{};
    MaterialProperties properties_{};
    BoundaryConditions conditions_{};
};
