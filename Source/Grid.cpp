#include "Grid.h"
#include "InputData.h"

Grid::Grid(const InputData& inpD, const NodesData& nD)
    : Nr_(nD.Nodes_r_)
    , Nz_(nD.Nodes_z_)
    , Ri_(inpD.geometry().innerRadius)
    , Ro_(inpD.geometry().outerRadius)
    , H_(inpD.geometry().height)
{
    dr_ = (Ro_ - Ri_) / (Nr_ - 1);
    dz_ = H_ / (Nz_ - 1);
}

std::size_t Grid::convertIndices(std::size_t i, std::size_t j) const
{
    return { Nr_ * j + i };
}
