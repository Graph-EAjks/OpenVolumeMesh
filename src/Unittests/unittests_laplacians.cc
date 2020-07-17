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





class DualLaplacianTest : public TetrahedralMeshBase{

public:

    //generate a basic mesh for further tests
    void SetUp() override{
        generateMinimalDiamondTetrahedralMesh(mesh_);
    }

private:

    void generateMinimalDiamondTetrahedralMesh(TetrahedralMesh& mesh){

        mesh.clear();

        double two_pi_5 = 2 * M_PI / 5.f;
        //pentagon shape
        Vec3d p0 = {cos(0 * two_pi_5) , sin(0 * two_pi_5), 0};
        Vec3d p1 = {cos(1 * two_pi_5) , sin(1 * two_pi_5), 0};
        Vec3d p2 = {cos(2 * two_pi_5) , sin(2 * two_pi_5), 0};
        Vec3d p3 = {cos(3 * two_pi_5) , sin(3 * two_pi_5), 0};
        Vec3d p4 = {cos(4 * two_pi_5) , sin(4 * two_pi_5), 0};

        Vec3d p5 = {0,0,1};
        Vec3d p6 = {0,0,0};

        VertexHandle v0 = mesh.add_vertex(p0);
        VertexHandle v1 = mesh.add_vertex(p1);
        VertexHandle v2 = mesh.add_vertex(p2);
        VertexHandle v3 = mesh.add_vertex(p3);
        VertexHandle v4 = mesh.add_vertex(p4);
        VertexHandle v5 = mesh.add_vertex(p5);
        VertexHandle v6 = mesh.add_vertex(p6);


        //'top faces
        FaceHandle f0 = mesh.add_face({v0, v1, v5});
        FaceHandle f1 = mesh.add_face({v1, v2, v5});
        FaceHandle f2 = mesh.add_face({v2, v3, v5});
        FaceHandle f3 = mesh.add_face({v3, v4, v5});
        FaceHandle f4 = mesh.add_face({v4, v0, v5});

        //bottom faces
        FaceHandle f5 = mesh.add_face({v0, v1, v6});
        FaceHandle f6 = mesh.add_face({v1, v2, v6});
        FaceHandle f7 = mesh.add_face({v2, v3, v6});
        FaceHandle f8 = mesh.add_face({v3, v4, v6});
        FaceHandle f9 = mesh.add_face({v4, v0, v6});

        //interior faces
        FaceHandle f10 = mesh.add_face({v0, v5, v6});
        FaceHandle f11 = mesh.add_face({v1, v5, v6});
        FaceHandle f12 = mesh.add_face({v2, v5, v6});
        FaceHandle f13 = mesh.add_face({v3, v5, v6});
        FaceHandle f14 = mesh.add_face({v4, v5, v6});

        //and cells
        mesh.add_cell({mesh.halfface_handle(f0, 0),
                       mesh.halfface_handle(f5, 1),
                       mesh.halfface_handle(f10,0),
                       mesh.halfface_handle(f11,1)});

        mesh.add_cell({mesh.halfface_handle(f1, 0),
                       mesh.halfface_handle(f6, 1),
                       mesh.halfface_handle(f11,0),
                       mesh.halfface_handle(f12,1)});

        mesh.add_cell({mesh.halfface_handle(f2, 0),
                       mesh.halfface_handle(f7, 1),
                       mesh.halfface_handle(f12,0),
                       mesh.halfface_handle(f13,1)});

        mesh.add_cell({mesh.halfface_handle(f3, 0),
                       mesh.halfface_handle(f8, 1),
                       mesh.halfface_handle(f13,0),
                       mesh.halfface_handle(f14,1)});

        mesh.add_cell({mesh.halfface_handle(f4, 0),
                       mesh.halfface_handle(f9, 1),
                       mesh.halfface_handle(f14,0),
                       mesh.halfface_handle(f10,1)});

        ASSERT_EQ(mesh.n_cells(), 5);

    }

};



TEST_F(DualLaplacianTest, CreateDualLaplacian){

    Laplacian<DualLaplacian, TetrahedralMesh> laplacian(mesh_);

}


TEST_F(DualLaplacianTest, GetPerHalfedgeWeight){

    Laplacian<DualLaplacian, TetrahedralMesh> laplacian(mesh_);

   ASSERT_DOUBLE_EQ(laplacian.halfedge_weight(HalfEdgeHandle(30)), 0.908178);
}


#if 0
TEST_F(DualLaplacianTest, GetPerVertexLaplacian){

    Laplacian<DualLaplacian, TetrahedralMesh> laplacian(mesh_);

    for(auto v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++){

        std::cout<<" laplacian for vertex "<<*v_it<<" : "<<laplacian[*v_it]<<std::endl;
    }

}

#endif

