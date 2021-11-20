#pragma once

#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <array>

namespace OpenVolumeMesh {

/// Compute area of a triangle face from vertex positions.
template<typename MeshT>
double compute_triangle_area(MeshT const&_mesh, FH _fh)
{
    using Point = typename MeshT::Point;
    std::array<Point, 3> pos;
    int idx = 0;
    for (auto vh: _mesh.face_vertices(_fh)) {
        if (idx == 3) {
            assert(false && "Face has more than 3 vertices!");
            break;
        }
        pos[idx++] = _mesh.vertex(vh);
    }
    const auto &[a,b,c] = pos;
    // TODO: check and improve numeric stability, maybe a variant of Heron's formula is better?
    return .5 * cross(b-a, c-a).norm();
}

} // namespace OpenVolumeMesh
