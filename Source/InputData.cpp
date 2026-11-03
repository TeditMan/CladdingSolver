#include "InputData.h"

#include <cmath>
#include <stdexcept>

InputData::InputData(Pressure pressure, Geometry geometry, MaterialProperties properties)
    : pressure_(pressure)
    , geometry_(geometry)
    , properties_(properties)
{
    if (!std::isfinite(pressure.Pgap) ||
        !std::isfinite(pressure.Pcool) ||
        !(pressure.Pgap > 0.0) ||
        !(pressure.Pcool > 0.0))
    {
        throw std::invalid_argument{"Invalid pressure values"};
    }

    if(!(std::isfinite(geometry.innerRadius)) ||
       !(std::isfinite(geometry.outerRadius)) ||
       !(std::isfinite(geometry.height)) ||
       !(geometry.innerRadius > 0.0) ||
       !(geometry.outerRadius > geometry.innerRadius) ||
       !(geometry.height > 0.0))
    {
        throw std::invalid_argument{"Invalid geometry"};
    }

    if (!std::isfinite(properties.youngModulus) ||
        !std::isfinite(properties.poissonRatio) ||
        !(properties.youngModulus > 0.0) ||
        !(properties.poissonRatio > -1.0 && properties.poissonRatio < 0.5))
    {
        throw std::invalid_argument{"Invalid material properties"};
    }

    conditions_.simga_rr_Ri = -pressure_.Pgap;
    conditions_.simga_rr_Ro = -pressure_.Pcool;
    conditions_.sigma_zz_H =
        (pressure.Pgap * geometry.innerRadius * geometry.innerRadius -
         pressure.Pcool * geometry.outerRadius * geometry.outerRadius) /
        (geometry.outerRadius * geometry.outerRadius -
         geometry.innerRadius * geometry.innerRadius);
}
