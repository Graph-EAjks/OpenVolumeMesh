#pragma once

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Attribs/TetHeightAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <array>

namespace OpenVolumeMesh {

/// compute integrated divergence of a per-cell 3-D vector field
/// on the 1-ball (set of incident cells) around a vertex.
template<typename MeshT, typename Vec>
typename MeshT::Point::value_type compute_vertex_divergence(
    MeshT const&_mesh,
    NormalAttrib<MeshT> const &_normal,
    TriangleAreaAttrib<MeshT> const &_area,
    PropertyPtr<Vec, Entity::Cell> const &_f,
    VH _vh)
{
    using Scalar = typename MeshT::Point::value_type;
    Scalar div = 0.;
    for (const auto ch: _mesh.vertex_cells(vh)) {
        TetTopology topo{_mesh, ch, vh};
        auto opp = topo.hfh<TetTopology::OppA>;
        div += dot(_f[ch], _normal[opp.opposite_handle()]) * _area[opp.face_handle()];
    }
    return div;
}
} // namespace OpenVolumeMesh
