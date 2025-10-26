#pragma once

#include <OpenVolumeMesh/Core/Iterators.hh>
#include <OpenVolumeMesh/Unstable/SmartHandles.hh>

namespace OpenVolumeMesh {
template<typename Derived, typename Parent>
class AddSmartIterators : public Parent
{
private:
    const Derived* derived() const { return static_cast<const Derived*>(this);}

protected:

    template <class Circulator>
    static Circulator make_end_circulator(const Circulator& _circ)
    {
        Circulator end = _circ;
        if (end.valid()) {
            end.lap(_circ.max_laps());
            end.valid(false);
        }
        return end;
    }


public:
    template<typename HandleT>
    auto make_smart(HandleT _h) const {
        return OpenVolumeMesh::make_smart(_h, *derived());
    }

    SmartRange<const Derived, VertexIter> vertices() const {
        return make_smart_range(derived(), derived()->vertices_begin(), derived()->vertices_end());
    }
    SmartRange<const Derived, EdgeIter> edges() const {
        return make_smart_range(derived(), derived()->edges_begin(), derived()->edges_end());
    }
    SmartRange<const Derived, HalfEdgeIter> halfedges() const {
        return make_smart_range(derived(), derived()->halfedges_begin(), derived()->halfedges_end());
    }
    SmartRange<const Derived, FaceIter> faces() const {
        return make_smart_range(derived(), derived()->faces_begin(), derived()->faces_end());
    }
    SmartRange<const Derived, HalfFaceIter> halffaces() const {
        return make_smart_range(derived(), derived()->halffaces_begin(), derived()->halffaces_end());
    }
    SmartRange<const Derived, CellIter> cells() const {
        return make_smart_range(static_cast<const Derived*>(this), derived()->cells_begin(), derived()->cells_end());
    }
    SmartRange<const Derived, VertexVertexIter> vertex_vertices(VertexHandle _h, int _max_laps = 1) const {
        VertexVertexIter begin = derived()->vv_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexIHalfEdgeIter> incoming_halfedges(VertexHandle _h, int _max_laps = 1) const {
        VertexIHalfEdgeIter begin = derived()->vih_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexOHalfEdgeIter> outgoing_halfedges(VertexHandle _h, int _max_laps = 1) const {
        VertexOHalfEdgeIter begin = derived()->voh_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexEdgeIter> vertex_edges(VertexHandle _h, int _max_laps = 1) const {
        VertexEdgeIter begin = derived()->ve_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexHalfFaceIter> vertex_halffaces(VertexHandle _h, int _max_laps = 1) const {
        VertexHalfFaceIter begin = derived()->vhf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexCellIter> vertex_cells(VertexHandle _h, int _max_laps = 1) const {
        VertexCellIter begin = derived()->vc_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfEdgeHalfFaceIter> halfedge_halffaces(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeHalfFaceIter begin = derived()->hehf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfEdgeFaceIter> halfedge_faces(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeFaceIter begin = derived()->hef_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfEdgeCellIter> halfedge_cells(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeCellIter begin = derived()->hec_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, EdgeHalfFaceIter> edge_halffaces(EdgeHandle _h, int _max_laps = 1) const {
        EdgeHalfFaceIter begin = derived()->ehf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, EdgeFaceIter> edge_faces(EdgeHandle _h, int _max_laps = 1) const {
        EdgeFaceIter begin = derived()->ef_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, EdgeCellIter> edge_cells(EdgeHandle _h, int _max_laps = 1) const {
        EdgeCellIter begin = derived()->ec_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfFaceHalfEdgeIter> halfface_halfedges(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceHalfEdgeIter begin = derived()->hfhe_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfFaceEdgeIter> halfface_edges(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceEdgeIter begin = derived()->hfe_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, FaceVertexIter> face_vertices(FaceHandle _h, int _max_laps = 1) const {
        FaceVertexIter begin = derived()->fv_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, FaceHalfEdgeIter> face_halfedges(FaceHandle _h, int _max_laps = 1) const {
        FaceHalfEdgeIter begin = derived()->fhe_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, FaceEdgeIter> face_edges(FaceHandle _h, int _max_laps = 1) const {
        FaceEdgeIter begin = derived()->fe_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellVertexIter> cell_vertices(CellHandle _h, int _max_laps = 1) const {
        CellVertexIter begin = derived()->cv_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellHalfEdgeIter> cell_halfedges(CellHandle _h, int _max_laps = 1) const {
        CellHalfEdgeIter begin = derived()->che_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellEdgeIter> cell_edges(CellHandle _h, int _max_laps = 1) const {
        CellEdgeIter begin = derived()->ce_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellHalfFaceIter> cell_halffaces(CellHandle _h, int _max_laps = 1) const {
        CellHalfFaceIter begin = derived()->chf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellFaceIter> cell_faces(CellHandle _h, int _max_laps = 1) const {
        CellFaceIter begin = derived()->cf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, CellCellIter> cell_cells(CellHandle _h, int _max_laps = 1) const {
        CellCellIter begin = derived()->cc_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, HalfFaceVertexIter> halfface_vertices(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceVertexIter begin = derived()->hfv_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, BoundaryHalfFaceHalfFaceIter> boundary_halfface_halffaces(HalfFaceHandle _h, int _max_laps = 1) const {
        BoundaryHalfFaceHalfFaceIter begin = derived()->bhfhf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
    SmartRange<const Derived, VertexFaceIter> vertex_faces(VertexHandle _h, int _max_laps = 1) const {
        VertexFaceIter begin = derived()->vf_iter(_h, _max_laps);
        return make_smart_range(derived(), begin, derived()->make_end_circulator(begin));
    }
};
} // namespace OpenVolumeMesh
