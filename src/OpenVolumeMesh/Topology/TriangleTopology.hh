#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>

namespace OpenVolumeMesh {

class OVM_EXPORT TriangleTopology
{
public:
    enum VertexLabel : uint8_t {
        A, B, C
    };
    enum HalfEdgeLabel : uint8_t {
        AB, BC, CA
    };
    TriangleTopology(TetrahedralMeshTopologyKernel const &mesh, FH fh);

    template<VertexLabel I>
    constexpr VH vh() const {return vh_[I];}

    constexpr VH a() const {return vh<A>();}
    constexpr VH b() const {return vh<B>();}
    constexpr VH c() const {return vh<C>();}

    template<HalfEdgeLabel I>
    constexpr HEH heh() const {return heh_[I];}

    constexpr HEH ab() const {return heh<AB>();}
    constexpr HEH bc() const {return heh<BC>();}
    constexpr HEH ca() const {return heh<CA>();}
private:
    std::array<VH, 3> vh_;
    std::array<HEH, 3> heh_;
};

} // namespace OpenVolumeMesh

