//
// Created by Denis Kalmykov on 29.07.22.
//
#include "unittests_common.hh"
#include "OpenVolumeMesh/FileManager/FileManager.hh"
#include "OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh"


class SplitIteratorTest : public ::testing::Test {
    using TetMesh = OpenVolumeMesh::TetrahedralGeometryKernel<OpenVolumeMesh::Vec3d, OpenVolumeMesh::TetrahedralMeshTopologyKernel>;
protected:
    void SetUp() override {
        inFileName = "/Users/dkalmykov/Development/OpenVolumeMesh/src/Unittests/TestFiles/5.ovm";
        fm.readFile(inFileName, tetmesh, true, true);
    }
    // void TearDown() override {}

    OpenVolumeMesh::IO::FileManager fm;
    TetMesh tetmesh;
    std::string inFileName;
};

TEST_F(SplitIteratorTest, IterationClassicsTest) {
    auto edge_ctr = 0;
    for (auto eh = tetmesh.edges_begin(); eh != tetmesh.edges_end(); ++eh) {
//            std::cout << eh->idx() << std::endl;
        edge_ctr++;
        auto edge_vertices = tetmesh.edge_vertices(*eh);
        if (tetmesh.is_boundary(edge_vertices[0]) and tetmesh.is_boundary(edge_vertices[1]) and
            !tetmesh.is_boundary(*eh)) {
            auto new_vertex = tetmesh.split_edge(*eh);
        }
    }
//        std::cout << "----------" << std::endl;
//        std::cout << edge_ctr << std::endl;
//        std::cout << "----------" << std::endl;
    tetmesh.collect_garbage();
    EXPECT_EQ(0, 0);
}

TEST_F(SplitIteratorTest, IterationForEachTest) {
    auto edge_ctr = 0;
    for(const auto &eh: tetmesh.edges()){
//        std::cout << eh.idx() << std::endl;
        edge_ctr++;
        auto edge_vertices = tetmesh.edge_vertices(eh);
        if(tetmesh.is_boundary(edge_vertices[0]) and tetmesh.is_boundary(edge_vertices[1]) and !tetmesh.is_boundary(eh)){
            auto new_vertex = tetmesh.split_edge(eh);
        }
    }
    std::cout << "----------" << std::endl;
    std::cout << edge_ctr << std::endl;
    std::cout << "----------" << std::endl;
    tetmesh.collect_garbage();
    EXPECT_EQ(0, 0);
}
