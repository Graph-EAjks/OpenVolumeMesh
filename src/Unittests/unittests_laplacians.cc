#include "unittests_common.hh"
#include <OpenVolumeMesh/Mesh/PolyhedralMeshLaplacians.hh>


using namespace OpenVolumeMesh;


class PolyhedralMeshLaplacianTest : public PolyhedralMeshBase{

public:

    //generate a basic mesh for further tests
    void SetUp() override{
        generatePolyhedralMesh(mesh_);
    }

private:

};

TEST_F(PolyhedralMeshLaplacianTest, CreateLaplacian){

    PolyhedralMeshLaplacian<PolyhedralMesh> laplacian(mesh_);

}



TEST_F(PolyhedralMeshLaplacianTest, EdgeUniformLaplacianEqualsOneEverywhere){

    PolyhedralMeshLaplacian<PolyhedralMesh> laplacian(mesh_);
    ASSERT_GT(mesh_.e_vertices(), 0);

    for(auto e_it = mesh_.edges_begin(); e_it != mesh_.edges_end(); e_it++){
        ASSERT_EQ(laplacian[*e_it], 1);

    }

}



TEST_F(PolyhedralMeshLaplacianTest, VertexUniformLaplacianEqualsAverageOfNeighborsPositions){

    PolyhedralMeshLaplacian<PolyhedralMesh> laplacian(mesh_);
    ASSERT_GT(mesh_.n_vertices(), 0);

    for(auto v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++){

        Vec3d average_position({0,0,0});
        double vertex_degree(0);
        auto vv_it = mesh_.vv_iter(*v_it);
        while(vv_it.valid()){
            average_position += mesh_.vertex(*vv_it);
            vertex_degree++;

            vv_it++;
        }
        average_position /= vertex_degree;


        ASSERT_EQ(laplacian[*v_it], average_position);

    }

}
