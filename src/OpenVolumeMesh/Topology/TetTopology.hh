#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
//#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


class OVM_EXPORT TetTopology {
public:
    enum VertexLabel : uint8_t {
        A, B, C, D
    };
    enum HalfEdgeLabel : uint8_t {
        // halfedge i+3 is vertex-disjunct with halfedge i (\in [0..2])
        AB, BC, CA,
        CD, AD, BD,
    };
    enum OppHalfEdgeLabel : uint8_t {
        BA = HalfEdgeLabel::AB,
        CB = HalfEdgeLabel::BC,
        AC = HalfEdgeLabel::CA,
        DC = HalfEdgeLabel::CD,
        DA = HalfEdgeLabel::AD,
        DB = HalfEdgeLabel::BD,
    };
    enum HalfFaceLabel : uint8_t {
        BDC = 0, DCB = 0, CBD = 0, OppA = 0,
        ACD = 1, CDA = 1, DAC = 1, OppB = 1,
        BAD = 2, ADB = 2, DBA = 2, OppC = 2,
        ABC = 3, BCA = 3, CAB = 3, OppD = 3,
    };

    /// Find a vertex labeling for the cell `ch` with the given
    /// halfface abc and vertex a.
    ///
    /// `abc` must be a a halfface of a the topological tetrahedron `ch`,
    /// and contain the vertex `a`.
    TetTopology(TopologyKernel const &mesh, CH ch, HFH abc, VH a=VH());

    /// Find a vertex labeling for the cell incident to `abc`, such that
    /// `a` is the given vertex handle, and a, b, and c form `abc`.
    ///
    /// `ch` MUST be a topological tetrahedron
    /// `hfh` MUST be a halface of a topological tetrahedron and contain `a`.
    TetTopology(TopologyKernel const &mesh, HFH abc, VH a=VH());

    /// Find a vertex labeling for the cell `ch`, such that `a` is the given
    /// vertex.
    /// `ch` MUST be a topological tetrahedron
    /// `a` MUST be a vertex of `ch`.
    TetTopology(TopologyKernel const &mesh, CH ch, VH a);

    /// Find any vertex labeling for the tetrahedral `ch`.
    /// `ch` MUST be a topological tetrahedron.
    TetTopology(TopologyKernel const &mesh, CH ch);

    bool operator==(TetTopology const& other) const;

    const auto &halfface_handles() const {return hfh_;}

    template<VertexLabel I>
    constexpr VH vh() const {return vh_[I];}
    constexpr VH vh(VertexLabel i) const {return vh_[i];}
    constexpr VH vh(int i) const {return vh_[i];}

    template<HalfEdgeLabel I>
    constexpr HEH heh() const {return heh_[I];}
    constexpr HEH heh(HalfEdgeLabel i) const {return heh_[i];}

    template<OppHalfEdgeLabel I>
    constexpr HEH heh() const {return heh_[I].opposite_handle();}
    constexpr HEH heh(OppHalfEdgeLabel i) const {return heh_[i].opposite_handle();}

    /// hfh<i>() is opposite to the vertex obtained by corner(i)
    template<HalfFaceLabel I>
    constexpr HFH hfh() const {return hfh_[I];}
    constexpr HFH hfh(HalfFaceLabel i) const {return hfh_[i];}
    constexpr HFH hfh(int i) const {return hfh_[i];}

    template<unsigned int X, unsigned int Y>
    constexpr HEH heh() const {
        static_assert(X<4);
        static_assert(Y<4);
        static_assert(X!=Y);
        if constexpr (X == 0 && Y == 1) return ab();
        if constexpr (X == 0 && Y == 2) return ac();
        if constexpr (X == 0 && Y == 3) return ad();

        if constexpr (X == 1 && Y == 0) return ba();
        if constexpr (X == 1 && Y == 2) return bc();
        if constexpr (X == 1 && Y == 3) return bd();

        if constexpr (X == 2 && Y == 0) return ca();
        if constexpr (X == 2 && Y == 1) return cb();
        if constexpr (X == 2 && Y == 3) return cd();

        if constexpr (X == 3 && Y == 0) return da();
        if constexpr (X == 3 && Y == 1) return db();
        if constexpr (X == 3 && Y == 2) return dc();
    }

    constexpr VH a() const {return vh<A>();}
    constexpr VH b() const {return vh<B>();}
    constexpr VH c() const {return vh<C>();}
    constexpr VH d() const {return vh<D>();}

    // half-edges
    constexpr HEH ab() const {return heh<AB>();}
    constexpr HEH ca() const {return heh<CA>();}
    constexpr HEH ad() const {return heh<AD>();}
    constexpr HEH bc() const {return heh<BC>();}
    constexpr HEH bd() const {return heh<BD>();}
    constexpr HEH cd() const {return heh<CD>();}

    constexpr HEH ac() const {return ca().opposite_handle();}
    constexpr HEH ba() const {return ab().opposite_handle();}
    constexpr HEH da() const {return ad().opposite_handle();}
    constexpr HEH dc() const {return cd().opposite_handle();}
    constexpr HEH db() const {return bd().opposite_handle();}
    constexpr HEH cb() const {return bc().opposite_handle();}

    // half-faces with all permutations
    constexpr HFH bdc() const {return hfh_[0];}
    constexpr HFH dcb() const {return hfh_[0];}
    constexpr HFH cbd() const {return hfh_[0];}

    constexpr HFH acd() const {return hfh_[1];}
    constexpr HFH cda() const {return hfh_[1];}
    constexpr HFH dac() const {return hfh_[1];}

    constexpr HFH bad() const {return hfh_[2];}
    constexpr HFH adb() const {return hfh_[2];}
    constexpr HFH dba() const {return hfh_[2];}

    constexpr HFH abc() const {return hfh_[3];}
    constexpr HFH bca() const {return hfh_[3];}
    constexpr HFH cab() const {return hfh_[3];}

private:
    std::array<VH, 4> vh_;
    std::array<HEH, 6> heh_;
    std::array<HFH, 4> hfh_;
};

} // namespace OpenVolumeMesh
