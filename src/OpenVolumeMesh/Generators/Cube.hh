#pragma once
#include <OpenVolumeMesh/Core/Handles.hh>

namespace OpenVolumeMesh {

template<typename MeshT>
FaceHandle add_axis_aligned_hex_cell(
        MeshT &_mesh,
        typename MeshT::Point const & _position,
        double _length)
{
    using Vector = typename MeshT::Point;
    int object_id;

    const double halfSize = 0.5*_length;
    std::vector<VertexHandle> vertices(8);
    // TODO: verify that the order matches the convention documented in HexahedralTopologyKernel
    vertices[0] = _mesh.add_vertex(Vector(-halfSize, -halfSize, halfSize)+_position);
    vertices[1] = _mesh.add_vertex(Vector( halfSize, -halfSize, halfSize)+_position);
    vertices[2] = _mesh.add_vertex(Vector( halfSize,  halfSize, halfSize)+_position);
    vertices[3] = _mesh.add_vertex(Vector(-halfSize,  halfSize, halfSize)+_position);
    vertices[4] = _mesh.add_vertex(Vector(-halfSize, -halfSize,-halfSize)+_position);
    vertices[5] = _mesh.add_vertex(Vector(-halfSize,  halfSize,-halfSize)+_position);
    vertices[6] = _mesh.add_vertex(Vector( halfSize,  halfSize,-halfSize)+_position);
    vertices[7] = _mesh.add_vertex(Vector( halfSize, -halfSize,-halfSize)+_position);
    // TODO: need add_hex that works on Polyhedral + Hexa TopologyKernel

    return _mesh.add_cell(vertices);
}
// TODO: add_axis_aligned_hex_grid
// TODO: versions with transformation matrix (or dx, dy, dz vectors, in lieu of a matrix class)

} // namespace OpenVolumeMesh 
