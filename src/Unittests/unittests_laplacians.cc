#include "unittests_common.hh"
#include <OpenVolumeMesh/Mesh/PolyhedralMeshLaplacians.hh>


using namespace OpenVolumeMesh;


class PolyhedralMeshLaplacianTest : public PolyhedralMeshBase{

public:

private:

};

TEST_F(PolyhedralMeshLaplacianTest, CreateLaplacian){

    PolyhedralMeshLaplacian<PolyhedralMesh> laplacian(mesh_);

}



TEST_F(PolyhedralMeshLaplacianTest, EdgeUniformLaplacianEqualsOneEverywhere){

    PolyhedralMeshLaplacian<PolyhedralMesh> laplacian(mesh_);

    for(auto e_it = mesh_.edges_begin(); e_it != mesh_.edges_end(); e_it++){
        ASSERT_EQ(laplacian[*e_it], 1);

    }

}
