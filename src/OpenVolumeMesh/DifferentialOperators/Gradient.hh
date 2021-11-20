#pragma once

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Attribs/TetHeightAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <array>

namespace OpenVolumeMesh {

/// Compute gradient of a scalar per-vertex function `_f` in a tet.
/// template parameter `MeshT` must be a geometric tetrahedral mesh
/// TetHeightAttrib and NormalAttrib must be available and computed.
template<typename MeshT, typename PropT>
typename MeshT::Point compute_tet_gradient(
    MeshT const&_mesh,
    NormalAttrib<MeshT> const &_normal,
    TetHeightAttrib<MeshT> const &_height,
    PropT const &_f,
    CH _ch
    )
{
    using Point = typename MeshT::Point;

    TetTopology topo{_mesh, _ch};

    std::array<HFH, 4> opp {
        topo.hfh<TetTopology::OppA>(),
        topo.hfh<TetTopology::OppB>(),
        topo.hfh<TetTopology::OppC>(),
        topo.hfh<TetTopology::OppD>()};

    Point grad {0., 0., 0.};
    for (int i = 0; i < 4; ++i) {
        grad += (_f[topo.vh(i)] / _height[opp[i]]) * _normal[opp[i]];
    }
    return grad;
}
} // namespace OpenVolumeMesh
