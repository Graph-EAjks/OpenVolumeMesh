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

    UniformLaplacian<PolyhedralMesh> laplacian(mesh_);

}



TEST_F(PolyhedralMeshLaplacianTest, HalfEdgeUniformLaplacianEqualsOneEverywhere){

    UniformLaplacian<PolyhedralMesh> laplacian(mesh_);
    ASSERT_GT(mesh_.n_halfedges(), 0);

    for(auto h_it = mesh_.halfedges_begin(); h_it != mesh_.halfedges_end(); h_it++){
        ASSERT_EQ(laplacian.halfedge_weight(*h_it), 1);

    }
}

TEST_F(PolyhedralMeshLaplacianTest, AccessHaldedgeWeightWithSubscriptOperator){

    Laplacian<UniformLaplacian, PolyhedralMesh> laplacian(mesh_);
    ASSERT_GT(mesh_.n_halfedges(), 0);

    for(auto h_it = mesh_.halfedges_begin(); h_it != mesh_.halfedges_end(); h_it++){
        ASSERT_EQ(laplacian[*h_it], 1);

    }

}


TEST_F(PolyhedralMeshLaplacianTest, VertexUniformLaplacianEqualsAverageOfNeighborsPositions){

    Laplacian<UniformLaplacian, PolyhedralMesh> laplacian(mesh_);
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



template<class _polyhedral_mesh>
class CustomLaplacian : BaseLaplacian<_polyhedral_mesh>{
public:
    CustomLaplacian(_polyhedral_mesh& mesh) : BaseLaplacian<_polyhedral_mesh>(mesh){}

    double halfedge_weight(const HalfEdgeHandle& edge) const { return 0;}
};



TEST_F(PolyhedralMeshLaplacianTest, CreateCustomLaplacians){

    Laplacian<CustomLaplacian, PolyhedralMesh> laplacian(mesh_);
}

