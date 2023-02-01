#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>

namespace OpenVolumeMesh {

template<class TopoKernel>
class SmartVertexHandle {
public:
    SmartVertexHandle(VH _vh, TopoKernel const&_kernel)
        : vh_(_vh)
        , kernel_(_kernel)
    {}
    operator VH(){ return vh_; }
    auto outgoing_halfedges() const { return kernel_.outgoing_halfedges(vh_); }
    auto incoming_halfedges() const { return kernel_.incoming_halfedges(vh_); }
    auto point() const { return kernel_.vertex(vh_); }
    auto faces() const { return kernel_.vertex_faces(vh_); }
    auto cells() const { return kernel_.vertex_cells(vh_); }
private:
    VH vh_;
    TopoKernel const& kernel_;
};

template<class TopoKernel>
class SmartEdgeHandle {
public:
    SmartEdgeHandle(EH _eh, TopoKernel const&_kernel)
        : eh_(_eh)
        , kernel_(_kernel)
    {}
    operator EH(){ return eh_; }
    auto faces() const { return kernel_.edge_faces(eh_); }
    auto cells() const { return kernel_.edge_cells(eh_); }
    auto h0() const { return make_smart(eh_.halfedge_handle(0), kernel_);}
    auto h1() const { return make_smart(eh_.halfedge_handle(1), kernel_);}
private:
    EH eh_;
    TopoKernel const& kernel_;
};

template<class TopoKernel>
class SmartHalfEdgeHandle {
public:
    SmartHalfEdgeHandle(HEH _heh, TopoKernel const&_kernel)
        : heh_(_heh)
        , kernel_(_kernel)
    {}
    operator HEH(){ return heh_; }
    auto from() const { return kernel_.from_vertex_handle(heh_); }
    auto to() const { return kernel_.to_vertex_handle(heh_); }
    auto edge() const { return make_smart(heh_.edge_handle(), kernel_);}
    auto faces() const { return kernel_.halfedge_faces(heh_); }
    auto halffaces() const { return kernel_.halfedge_halffaces(heh_); }
    auto cells() const { return kernel_.halfedge_cells(heh_); }
private:
    HEH heh_;
    TopoKernel const& kernel_;
};

template<class TopoKernel>
class SmartFaceHandle {
public:
    SmartFaceHandle(FH _fh, TopoKernel const&_kernel)
        : fh_(_fh)
        , kernel_(_kernel)
    {}
    operator FH(){ return fh_; }
    auto h0() const { return make_smart(fh_.halfface_handle(0), kernel_);}
    auto h1() const { return make_smart(fh_.halfface_handle(1), kernel_);}
    auto vertices() const { return kernel_.face_vertices(fh_); }
private:
    FH fh_;
    TopoKernel const& kernel_;
};

template<class TopoKernel>
class SmartHalfFaceHandle {
public:
    SmartHalfFaceHandle(HFH _hfh, TopoKernel const&_kernel)
        : hfh_(_hfh)
        , kernel_(_kernel)
    {}
    operator HFH(){ return hfh_; }
    auto from() const { return kernel_.from_vertex_handle(hfh_); }
    auto to() const { return kernel_.to_vertex_handle(hfh_); }
    auto face() const { return make_smart(hfh_.face_handle(), kernel_);}
    auto halfedges() const { return kernel_.halfedges(hfh_); }
    auto cell() const { return kernel_.incident_cell(hfh_); }
    auto vertices() const { return kernel_.halfface_vertices(hfh_); }
private:
    HFH hfh_;
    TopoKernel const& kernel_;
};

template<class TopoKernel>
class SmartCellHandle {
public:
    SmartCellHandle(CH _ch, TopoKernel const&_kernel)
        : ch_(_ch)
        , kernel_(_kernel)
    {}
    operator CH(){ return ch_; }
    auto halffaces() const { return kernel_.cell_halffaces(ch_); }
    auto vertices() const { return kernel_.cell_vertices(ch_); }
private:
    CH ch_;
    TopoKernel const& kernel_;
};

#if 0
template<typename Entity, class TopoKernel>
struct smart_handle {};

template<typename Entity, class TopoKernel>
using smart_handle_t = typename smart_handle<Entity,TopoKernel>::type;

template<class TopoKernel>
struct smart_handle<Entity::Cell, TopoKernel> { using type = SmartCellHandle<TopoKernel>;};
#endif

template<class TopoKernel>
auto make_smart(VH _h, TopoKernel const&_kernel) {
    return SmartVertexHandle{_h, _kernel};
}
template<class TopoKernel>
auto make_smart(EH _h, TopoKernel const&_kernel) {
    return SmartEdgeHandle{_h, _kernel};
}
template<class TopoKernel>
auto make_smart(HEH _h, TopoKernel const&_kernel) {
    return SmartHalfEdgeHandle{_h, _kernel};
}
template<class TopoKernel>
auto make_smart(FH _h, TopoKernel const&_kernel) {
    return SmartFaceHandle{_h, _kernel};
}
template<class TopoKernel>
auto make_smart(HFH _h, TopoKernel const&_kernel) {
    return SmartHalfFaceHandle{_h, _kernel};
}
template<class TopoKernel>
auto make_smart(CH _h, TopoKernel const&_kernel) {
    return SmartCellHandle{_h, _kernel};
}

} // namespace OpenVolumeMesh
