#include "unittests_common.hh"
#include <OpenVolumeMesh/Attribs/TetHeightAttrib.hh>
#include <OpenVolumeMesh/Attribs/TetVolumeAttrib.hh>
#include <OpenVolumeMesh/Attribs/TriangleAreaAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <OpenVolumeMesh/Geometry/Volume.hh>
#include <OpenVolumeMesh/DifferentialOperators/Gradient.hh>
#include <OpenVolumeMesh/Geometry/Area.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

namespace OpenVolumeMesh {

struct TetGeometryTest : public testing::Test {
  GeometricTetrahedralMeshV3d mesh;
  VH v0 = mesh.add_vertex({0, 0, 0});
  VH v1 = mesh.add_vertex({1, 0, 0});
  VH v2 = mesh.add_vertex({0, 1, 0});
  VH v3 = mesh.add_vertex({0, 0, 1});
  HFH hf0 = mesh.add_halfface(v0, v1, v2);
  HFH hf1 = mesh.add_halfface(v1, v0, v3);
  HFH hf2 = mesh.add_halfface(v2, v1, v3);
  HFH hf3 = mesh.add_halfface(v0, v2, v3);
  CH ch = mesh.add_cell({hf0, hf1, hf2, hf3});
  // TODO: apply random rotation & scaling and test co/contra/in-variance
};

TEST_F(TetGeometryTest, TetVolume){
    auto vol = compute_tet_volume(mesh, ch);
    EXPECT_DOUBLE_EQ(vol, 1./6);
    TetVolumeAttrib tva{mesh};
    tva.update();
    EXPECT_DOUBLE_EQ(tva[ch], vol);
}
TEST_F(TetGeometryTest, TriangleArea){
    auto a0 = compute_triangle_area(mesh, hf0.face_handle());
    EXPECT_DOUBLE_EQ(a0, .5);
    auto a1 = compute_triangle_area(mesh, hf1.face_handle());
    EXPECT_DOUBLE_EQ(a1, .5);
    auto a2 = compute_triangle_area(mesh, hf2.face_handle());
    EXPECT_DOUBLE_EQ(a2, sqrt(3./4));
    auto a3 = compute_triangle_area(mesh, hf3.face_handle());
    EXPECT_DOUBLE_EQ(a3, .5);
    TriangleAreaAttrib taa{mesh};
    taa.update();
    EXPECT_DOUBLE_EQ(taa[hf0.face_handle()], a0);
    EXPECT_DOUBLE_EQ(taa[hf1.face_handle()], a1);
    EXPECT_DOUBLE_EQ(taa[hf2.face_handle()], a2);
    EXPECT_DOUBLE_EQ(taa[hf3.face_handle()], a3);
}
TEST_F(TetGeometryTest, TetHeight){
    TetVolumeAttrib volume{mesh};
    TriangleAreaAttrib area{mesh};
    TetHeightAttrib height{mesh};
    volume.update();
    area.update();
    height.update(volume, area);
    EXPECT_TRUE(std::isnan(height[hf0.opposite_handle()]));
    EXPECT_TRUE(std::isnan(height[hf1.opposite_handle()]));
    EXPECT_TRUE(std::isnan(height[hf2.opposite_handle()]));
    EXPECT_TRUE(std::isnan(height[hf3.opposite_handle()]));
    EXPECT_DOUBLE_EQ(height[hf0], 1.);
    EXPECT_DOUBLE_EQ(height[hf1], 1.);
    EXPECT_DOUBLE_EQ(height[hf3], 1.);
    EXPECT_DOUBLE_EQ(height[hf2], sqrt(1./3)); // TODO
}

TEST_F(TetGeometryTest, TetGradient){

    NormalAttrib normal{mesh};
    normal.update_face_normals();
    TetVolumeAttrib volume{mesh};
    volume.update();
    TriangleAreaAttrib area{mesh};
    area.update();
    TetHeightAttrib height{mesh};
    height.update(volume, area);

    auto f = mesh.create_private_property<double, Entity::Vertex>();
    auto grad = [&](std::array<double, 4> _f) mutable {
        f[v0] = _f[0]; f[v1] = _f[1]; f[v2] = _f[2]; f[v3] = _f[3];
        return compute_tet_gradient(mesh, normal, height, f, ch);
    };

    Vec3d g;
    auto test_grad = [&]() {
        EXPECT_DOUBLE_EQ(f[v3] - f[v0], dot(g, mesh.vertex(v3) - mesh.vertex(v0)));
        EXPECT_DOUBLE_EQ(f[v2] - f[v0], dot(g, mesh.vertex(v2) - mesh.vertex(v0)));
        EXPECT_DOUBLE_EQ(f[v1] - f[v0], dot(g, mesh.vertex(v1) - mesh.vertex(v0)));

        EXPECT_DOUBLE_EQ(f[v3] - f[v1], dot(g, mesh.vertex(v3) - mesh.vertex(v1)));
        EXPECT_DOUBLE_EQ(f[v2] - f[v1], dot(g, mesh.vertex(v2) - mesh.vertex(v1)));

        EXPECT_DOUBLE_EQ(f[v3] - f[v2], dot(g, mesh.vertex(v3) - mesh.vertex(v2)));
    };
    g = grad({0, 0, 0, 1});
    EXPECT_DOUBLE_EQ(g[0], 0.);
    EXPECT_DOUBLE_EQ(g[1], 0.);
    EXPECT_DOUBLE_EQ(g[2], 1.);
    test_grad();

    g = grad({0, 0, 1, 0});
    EXPECT_DOUBLE_EQ(g[0], 0.);
    EXPECT_DOUBLE_EQ(g[1], 1.);
    EXPECT_DOUBLE_EQ(g[2], 0.);
    test_grad();

    g = grad({0, 1, 0, 0});
    EXPECT_DOUBLE_EQ(g[0], 1.);
    EXPECT_DOUBLE_EQ(g[1], 0.);
    EXPECT_DOUBLE_EQ(g[2], 0.);
    test_grad();

    g = grad({-1, 0, 0, 0});
    EXPECT_DOUBLE_EQ(g[0], 1.);
    EXPECT_DOUBLE_EQ(g[1], 1.);
    EXPECT_DOUBLE_EQ(g[2], 1.);
    test_grad();
}

} //namespace OpenVolumeMesh

