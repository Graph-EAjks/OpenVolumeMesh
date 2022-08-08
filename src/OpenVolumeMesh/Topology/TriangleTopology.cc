#include <OpenVolumeMesh/Topology/TriangleTopology.hh>

namespace OpenVolumeMesh {
TriangleTopology::TriangleTopology(
        const TetrahedralMeshTopologyKernel &_mesh,
        FH _fh)
{
    int idx = 0;
    for (const auto heh: _mesh.face_halfedges(_fh))
    {
        heh_[idx] = heh;
        vh_[idx] = _mesh.from_vertex_handle(heh);
        ++idx;
    }
    assert(idx == 3);
}

bool TriangleTopology::operator==(TriangleTopology const& other) const
{
    return vh_ == other.vh_
       && heh_ == other.heh_;
}

} // namespace OpenVolumeMesh
