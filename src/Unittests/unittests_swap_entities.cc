#include <iostream>

#include <gtest/gtest.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>

using namespace OpenVolumeMesh;

class SwapEntitiesTest : public testing::Test {
};

TEST_F(SwapEntitiesTest, DenisBug)
{
    Reader::VtkColorReader vtk_reader;
    GeometricTetrahedralMeshV3d mesh;
    bool success = vtk_reader.readFile("s04u_tetrahedron.vtk", mesh);
    ASSERT_EQ(success, true);

    auto ch = CH(27);

    const auto vhs = mesh.get_cell_vertices(ch);
    ASSERT_EQ(vhs.size(), 4);
    mesh.swap_vertex_indices(vhs[2],vhs[3]);
    ASSERT_EQ(vhs.size(), 4);
    const auto vhs2 = mesh.get_cell_vertices(ch);
    ASSERT_EQ(vhs2.size(), 4);
}

