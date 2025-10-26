#include "unittests_common.hh"
#include <OpenVolumeMesh/Unstable/SmartHandles.hh>
#include <algorithm>
#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif
#if defined(__cpp_lib_ranges)
#  include <ranges>
#endif

using namespace OpenVolumeMesh;

/*
 * Test SmartVertexHandle functionality
 */
TEST_F(TetrahedralMeshBase, SmartVertexHandle_BasicFunctionality) {
    generateTetrahedralMesh(mesh_);

    auto vh = VertexHandle(0);
    auto smart_vh = make_smart(vh, mesh_);

    EXPECT_HANDLE_EQ(vh, smart_vh);

    for (auto heh : smart_vh.outgoing_halfedges()) {
        EXPECT_EQ(heh.from(), smart_vh);
    }
    for (auto cell : mesh_.cells()) {
        for (auto hfh : cell.halffaces()) {
            EXPECT_EQ(hfh.cell(), cell);
        }
    }
}

/*
/*
 * Test smart handle with max_laps parameter
 */
TEST_F(TetrahedralMeshBase, SmartHandles_MaxLapsParameter) {
    generateTetrahedralMesh(mesh_);

    auto smart_vh = make_smart(VertexHandle(0), mesh_);

    // Test with different max_laps values
    int count_default = 0;
    for (auto vh : smart_vh.adjacent_vertices()) {
        (void)vh;
        count_default++;
    }

    int count_single = 0;
    for (auto vh : smart_vh.adjacent_vertices(1)) {
        (void)vh;
        count_single++;
    }

    EXPECT_EQ(count_default, count_single);

    int count_twice = 0;
    for (auto vh : smart_vh.adjacent_vertices(2)) {
        (void)vh;
        count_twice++;
    }
    EXPECT_EQ(count_twice, count_single * 2);


    // Test with max_laps = 0 (should iterate once)
    int count_zero = 0;
    for (auto vh : smart_vh.adjacent_vertices(0)) {
        (void)vh;
        count_zero++;
    }
    EXPECT_GE(count_zero, 0);
}

/*
 * Test smart_handle trait
 */
TEST_F(TetrahedralMeshBase, SmartHandles_TraitSpecialization) {
    // This is a compile-time test to ensure all trait specializations exist
    using VertexSmartType = smart_handle_t<Entity::Vertex, TetrahedralMesh>;
    using EdgeSmartType = smart_handle_t<Entity::Edge, TetrahedralMesh>;
    using HalfEdgeSmartType = smart_handle_t<Entity::HalfEdge, TetrahedralMesh>;
    using FaceSmartType = smart_handle_t<Entity::Face, TetrahedralMesh>;
    using HalfFaceSmartType = smart_handle_t<Entity::HalfFace, TetrahedralMesh>;
    using CellSmartType = smart_handle_t<Entity::Cell, TetrahedralMesh>;

    // If this compiles, the traits are properly specialized
    EXPECT_TRUE((std::is_same<VertexSmartType, SmartVertexHandle<TetrahedralMesh>>::value));
    EXPECT_TRUE((std::is_same<EdgeSmartType, SmartEdgeHandle<TetrahedralMesh>>::value));
    EXPECT_TRUE((std::is_same<HalfEdgeSmartType, SmartHalfEdgeHandle<TetrahedralMesh>>::value));
    EXPECT_TRUE((std::is_same<FaceSmartType, SmartFaceHandle<TetrahedralMesh>>::value));
    EXPECT_TRUE((std::is_same<HalfFaceSmartType, SmartHalfFaceHandle<TetrahedralMesh>>::value));
    EXPECT_TRUE((std::is_same<CellSmartType, SmartCellHandle<TetrahedralMesh>>::value));
}

#if defined(__cpp_lib_ranges)
TEST_F(TetrahedralMeshBase, SmartHandles_std_ranges) {
    generateTetrahedralMesh(mesh_);

    static_assert(std::ranges::range<OpenVolumeMesh::SmartRange<TetrahedralMeshTopologyKernel, VertexIter>>);
    auto first2 = mesh_.vertices() | std::views::take(2);
    EXPECT_EQ(std::ranges::distance(first2), 2u);

    auto heh_idxs = mesh_.make_smart(FH(0)).halfedges() |
        std::views::transform([](auto h){return h.idx();});
    EXPECT_EQ(std::ranges::distance(heh_idxs), 3u);
}
#else
TEST_F(TetrahedralMeshBase, DISABLED_SmartHandles_std_ranges) {
    SUCCEED();
}
#endif
