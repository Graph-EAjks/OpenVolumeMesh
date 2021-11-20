#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <array>

namespace OpenVolumeMesh {

/// Compute volume of a tet cell from the 3-D vertex embeddings.
/// template parameter `MeshT` must be a geometric tetrahedral mesh
template<typename MeshT>
double compute_tet_volume(MeshT const&_mesh, CH _ch)
{
    using Point = typename MeshT::Point;
    std::array<Point, 4> pos;
    int idx = 0;
    for (auto vh: _mesh.tet_vertices(_ch)) {
        pos[idx++] = _mesh.vertex(vh);
    }
    const auto &[a,b,c,d] = pos;
    return dot(b-a, cross(c-a, d-a)) * (1./6);
}

} // namespace OpenVolumeMesh
