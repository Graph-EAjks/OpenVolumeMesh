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

    auto testCH = CH(10783);

    std::cerr << "nc " << mesh.n_cells() << std::endl;

    const auto vhs = mesh.get_cell_vertices(testCH);
    for (const auto & vh : vhs) {
        std::cerr << vh << std::endl;
    }
    mesh.swap_vertex_indices(vhs[2],vhs[3]);
    const auto vhs2 = mesh.get_cell_vertices(testCH);
    for (const auto & vh : vhs2) {
        std::cerr << vh << std::endl;
    }
    std::cout << "CHECK!" << std::endl;

}

