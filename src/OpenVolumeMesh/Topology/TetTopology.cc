#include <OpenVolumeMesh/Topology/TetTopology.hh>

namespace OpenVolumeMesh {

static HFH find_halfface(TetrahedralMeshTopologyKernel const &_mesh, CH _ch, VH _vh) {
    for (const auto hfh: _mesh.cell_halffaces(_ch)) {
        for (const auto vh: _mesh.halfface_vertices(hfh)) {
            if (vh == _vh)
                return hfh;
        }
    }
    assert(false);
    return {};
}


TetTopology::TetTopology(const TetrahedralMeshTopologyKernel &mesh,
                         CH _ch,
                         HFH _abc,
                         VH _a)
{
    auto abc_he_it = mesh.hfhe_iter(_abc, 2);
    if (_a.is_valid()) {
        while(mesh.from_vertex_handle(*abc_he_it) != _a) {
            ++abc_he_it;
        }
    }
    heh_[AB] = *abc_he_it++;
    heh_[BC] = *abc_he_it++;
    heh_[CA] = *abc_he_it;

    vh_[A] = mesh.from_vertex_handle(ab());
    vh_[B] = mesh.from_vertex_handle(bc());
    vh_[C] = mesh.from_vertex_handle(ca());

    assert(b() == mesh.to_vertex_handle(ab()));
    assert(c() == mesh.to_vertex_handle(bc()));
    assert(a() == mesh.to_vertex_handle(ca()));

    hfh_[ABC] = _abc;
    for (const auto hfh: mesh.cell_halffaces(_ch)) {
        if (hfh == _abc) continue;
        auto heh_it = mesh.hfhe_iter(hfh, 2);
        while (heh_it) {
            auto heh = *(heh_it++);
            if (heh == ba()) {
                hfh_[BAD] = hfh;
                heh_[AD] = *heh_it;
                break;
            }
            if (heh == cb()) {
                hfh_[CBD] = hfh;
                heh_[BD] = *heh_it;
                break;
            }
            if (heh == ac()) {
                hfh_[ACD] = hfh;
                heh_[CD] = *heh_it;
                break;
            }
        }
    }
    vh_[D] = mesh.to_vertex_handle(ad());

    assert(d() != a());
    assert(d() != b());
    assert(d() != c());

    assert(d() == mesh.to_vertex_handle(ad()));
    assert(d() == mesh.to_vertex_handle(bd()));
    assert(d() == mesh.to_vertex_handle(cd()));

    assert(a() == mesh.from_vertex_handle(ad()));
    assert(b() == mesh.from_vertex_handle(bd()));
    assert(c() == mesh.from_vertex_handle(cd()));
}


TetTopology::TetTopology(TetrahedralMeshTopologyKernel const &_mesh,
        HFH _abc, VH _a)
    : TetTopology(_mesh, _mesh.incident_cell(_abc), _abc, _a)
{}

TetTopology::TetTopology(const TetrahedralMeshTopologyKernel &_mesh,
                         CH _ch, VH _a)
    : TetTopology(_mesh, _ch, find_halfface(_mesh, _ch, _a), _a)
{}

TetTopology::TetTopology(const TetrahedralMeshTopologyKernel &_mesh, CH _ch)
    : TetTopology(_mesh, *_mesh.chf_iter(_ch))
{}

bool TetTopology::operator==(TetTopology const& other) const
{
  return vh_ == other.vh_
    && heh_ == other.heh_
    && hfh_ == other.hfh_;
}


} // namespace OpenVolumeMesh
