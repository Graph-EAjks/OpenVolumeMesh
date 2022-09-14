/*
 * unittests_common.cc
 *
 *  Created on: Jan 10, 2012
 *      Author: kremer
 */

#include "unittests_common.hh"

using namespace OpenVolumeMesh::Geometry;

void PolyhedralMeshBase::generatePolyhedralMesh(PolyhedralMesh& _mesh) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    Vec3d p9(0.0, 0.0, 2.0);
    Vec3d p10(1.0, 0.0, 2.0);
    Vec3d p11(1.0, 1.0, 2.0);
    Vec3d p12(0.0, 1.0, 2.0);

    VertexHandle v1 = _mesh.add_vertex(p1);
    VertexHandle v2 = _mesh.add_vertex(p2);
    VertexHandle v3 = _mesh.add_vertex(p3);
    VertexHandle v4 = _mesh.add_vertex(p4);

    VertexHandle v5 = _mesh.add_vertex(p5);
    VertexHandle v6 = _mesh.add_vertex(p6);
    VertexHandle v7 = _mesh.add_vertex(p7);
    VertexHandle v8 = _mesh.add_vertex(p8);

    VertexHandle v9  = _mesh.add_vertex(p9);
    VertexHandle v10 = _mesh.add_vertex(p10);
    VertexHandle v11 = _mesh.add_vertex(p11);
    VertexHandle v12 = _mesh.add_vertex(p12);

    std::vector<VertexHandle> vertices;
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v3); vertices.push_back(v4);
    FaceHandle f1 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f2 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f3 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v4);
    vertices.push_back(v8); vertices.push_back(v5);
    FaceHandle f4 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v4); vertices.push_back(v3);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f5 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v6); vertices.push_back(v5);
    vertices.push_back(v1); vertices.push_back(v2);
    FaceHandle f6 = _mesh.add_face(vertices);

    // Add first cell
    std::vector<HalfFaceHandle> halffaces;
    halffaces.push_back(_mesh.halfface_handle(f1, 1)); halffaces.push_back(_mesh.halfface_handle(f2, 0));
    halffaces.push_back(_mesh.halfface_handle(f3, 1)); halffaces.push_back(_mesh.halfface_handle(f4, 1));
    halffaces.push_back(_mesh.halfface_handle(f5, 1)); halffaces.push_back(_mesh.halfface_handle(f6, 0));
    _mesh.add_cell(halffaces);

    vertices.clear();
    vertices.push_back(v9); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f7 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v10); vertices.push_back(v11);
    vertices.push_back(v7); vertices.push_back(v6);
    FaceHandle f8 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v8);
    vertices.push_back(v12); vertices.push_back(v9);
    FaceHandle f9 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v8); vertices.push_back(v7);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f10 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v10); vertices.push_back(v9);
    vertices.push_back(v5); vertices.push_back(v6);
    FaceHandle f11 = _mesh.add_face(vertices);

    halffaces.clear();
    halffaces.push_back(_mesh.halfface_handle(f2,  1)); halffaces.push_back(_mesh.halfface_handle(f7,  0));
    halffaces.push_back(_mesh.halfface_handle(f8,  1)); halffaces.push_back(_mesh.halfface_handle(f9,  1));
    halffaces.push_back(_mesh.halfface_handle(f10, 1)); halffaces.push_back(_mesh.halfface_handle(f11, 0));
    _mesh.add_cell(halffaces);
}

void HexahedralMeshBase::generateHexahedralMesh(HexahedralMesh& _mesh) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    Vec3d p9(0.0, 0.0, 2.0);
    Vec3d p10(1.0, 0.0, 2.0);
    Vec3d p11(1.0, 1.0, 2.0);
    Vec3d p12(0.0, 1.0, 2.0);

    VertexHandle v1 = _mesh.add_vertex(p1);
    VertexHandle v2 = _mesh.add_vertex(p2);
    VertexHandle v3 = _mesh.add_vertex(p3);
    VertexHandle v4 = _mesh.add_vertex(p4);

    VertexHandle v5 = _mesh.add_vertex(p5);
    VertexHandle v6 = _mesh.add_vertex(p6);
    VertexHandle v7 = _mesh.add_vertex(p7);
    VertexHandle v8 = _mesh.add_vertex(p8);

    VertexHandle v9  = _mesh.add_vertex(p9);
    VertexHandle v10 = _mesh.add_vertex(p10);
    VertexHandle v11 = _mesh.add_vertex(p11);
    VertexHandle v12 = _mesh.add_vertex(p12);

    std::vector<VertexHandle> vertices;
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v3); vertices.push_back(v4);
    FaceHandle f0 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f1 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f2 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v5);
    vertices.push_back(v8); vertices.push_back(v4);
    FaceHandle f3 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);
    vertices.push_back(v6); vertices.push_back(v5);
    FaceHandle f4 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v4); vertices.push_back(v3);
    vertices.push_back(v7); vertices.push_back(v8);
    FaceHandle f5 = _mesh.add_face(vertices);

    // Add first cell
    std::vector<HalfFaceHandle> halffaces;
    halffaces.push_back(_mesh.halfface_handle(f0, 1)); halffaces.push_back(_mesh.halfface_handle(f1, 0));
    halffaces.push_back(_mesh.halfface_handle(f2, 1)); halffaces.push_back(_mesh.halfface_handle(f3, 0));
    halffaces.push_back(_mesh.halfface_handle(f4, 0)); halffaces.push_back(_mesh.halfface_handle(f5, 1));
    _mesh.add_cell(halffaces);

    vertices.clear();
    vertices.push_back(v9); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f6 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v6); vertices.push_back(v10);
    vertices.push_back(v11); vertices.push_back(v7);
    FaceHandle f7 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v9);
    vertices.push_back(v12); vertices.push_back(v8);
    FaceHandle f8 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v5); vertices.push_back(v6);
    vertices.push_back(v10); vertices.push_back(v9);
    FaceHandle f9 = _mesh.add_face(vertices);

    vertices.clear();
    vertices.push_back(v8); vertices.push_back(v7);
    vertices.push_back(v11); vertices.push_back(v12);
    FaceHandle f10 = _mesh.add_face(vertices);

    halffaces.clear();
    halffaces.push_back(_mesh.halfface_handle(f1, 1)); halffaces.push_back(_mesh.halfface_handle(f6, 0));
    halffaces.push_back(_mesh.halfface_handle(f7, 1)); halffaces.push_back(_mesh.halfface_handle(f8, 0));
    halffaces.push_back(_mesh.halfface_handle(f9, 0)); halffaces.push_back(_mesh.halfface_handle(f10, 1));
    _mesh.add_cell(halffaces);
}

void TetrahedralMeshBase::generate_tetrahedral_mesh(TetrahedralMesh& _mesh) {

    _mesh.clear();

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    VertexHandle v1 = _mesh.add_vertex(p1);
    VertexHandle v2 = _mesh.add_vertex(p2);
    VertexHandle v3 = _mesh.add_vertex(p3);
    VertexHandle v4 = _mesh.add_vertex(p4);

    // Add  cell
    _mesh.add_cell(v1, v2, v3, v4);
}

void TetrahedralMeshBase::generate_tetrahedral_mesh_2(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }
    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3], true);
    _mesh.add_cell(vertices[0],  vertices[2], vertices[1], vertices[4]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_2T1V(TetrahedralMesh& _mesh) {

    _mesh.clear();
    VertexHandle vertices[7];
    for (int i = 0; i < 7; ++i) {
        vertices[i] = _mesh.add_vertex();
    }
    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    _mesh.add_cell(vertices[0], vertices[4], vertices[5], vertices[6]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_2T1E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[6];
    for (int i = 0; i < 6; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    _mesh.add_cell(vertices[0], vertices[1], vertices[4], vertices[5]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_6T1E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[8];
    for (int i = 0; i < 8; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    _mesh.add_cell(vertices[1], vertices[0], vertices[2], vertices[4]);
    _mesh.add_cell(vertices[0], vertices[1], vertices[5], vertices[6]);
    _mesh.add_cell(vertices[1], vertices[0], vertices[5], vertices[7]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_1T1F1E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3], true);
    std::vector<VertexHandle> face_vertices = {vertices[0], vertices[1], vertices[4]};
    _mesh.add_face(face_vertices);
}

void TetrahedralMeshBase::generate_non_manifold_tet_2F1E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 4; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<VertexHandle> faceVertices = {vertices[0], vertices[1], vertices[2]};
    _mesh.add_face(faceVertices);
    faceVertices = {vertices[0], vertices[1], vertices[3]};
    _mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generate_non_manifold_tet_3T1V3E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[7];
    for (int i = 0; i < 7; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[4]);
    _mesh.add_cell(vertices[0], vertices[3], vertices[1], vertices[5]);
    _mesh.add_cell(vertices[0], vertices[2], vertices[3], vertices[6]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_1V(TetrahedralMesh &_mesh) {
    _mesh.clear();

    _mesh.add_vertex();
}

void TetrahedralMeshBase::generate_non_manifold_tet_1E(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[2];
    for (int i = 0; i < 2; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_edge(vertices[0], vertices[1]);
}

void TetrahedralMeshBase::generate_non_manifold_tet_1F(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[3];
    for (int i = 0; i < 3; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<VertexHandle> faceVertices = {vertices[0], vertices[1], vertices[2]};
    _mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generate_non_manifold_tet_1T1F1V(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[6];
    for (int i = 0; i < 6; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    std::vector<VertexHandle> faceVertices = { vertices[0], vertices[4], vertices[5]};
    _mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generate_non_manifold_tet_2F1V(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<VertexHandle> vertices_face_1 = {vertices[0], vertices[1], vertices[2]};
    std::vector<VertexHandle> vertices_face_2 = {vertices[0], vertices[3], vertices[4]};
    _mesh.add_face(vertices_face_1);
    _mesh.add_face(vertices_face_2);
}

void TetrahedralMeshBase::generate_non_manifold_tet_1F1E1V(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 4; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<VertexHandle> face_vertices = {vertices[0], vertices[1], vertices[2]};
    _mesh.add_face(face_vertices);
    _mesh.add_edge(vertices[0], vertices[3]);
}

void TetrahedralMeshBase::generate_tet_without_cell(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 4; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_face({vertices[0], vertices[1], vertices[2]});
    _mesh.add_face({vertices[0], vertices[3], vertices[1]});
    _mesh.add_face({vertices[0], vertices[2], vertices[3]});
    _mesh.add_face({vertices[1], vertices[3], vertices[2]});
}

void TetrahedralMeshBase::generate_tri_without_face(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[3];
    for (int i = 0; i < 3; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_edge(vertices[0], vertices[1]);
    _mesh.add_edge(vertices[1], vertices[2]);
    _mesh.add_edge(vertices[2], vertices[0]);
}

void TetrahedralMeshBase::generate_tet_without_cells_and_faces(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_edge(vertices[0], vertices[1]);
    _mesh.add_edge(vertices[0], vertices[2]);
    _mesh.add_edge(vertices[0], vertices[3]);
    _mesh.add_edge(vertices[1], vertices[2]);
    _mesh.add_edge(vertices[1], vertices[3]);
    _mesh.add_edge(vertices[2], vertices[3]);
    _mesh.add_edge(vertices[1], vertices[4]);
    _mesh.add_edge(vertices[2], vertices[4]);
    _mesh.add_edge(vertices[3], vertices[4]);
}

void TetrahedralMeshBase::generate_tri_tet(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_edge(vertices[0], vertices[1]);
    _mesh.add_edge(vertices[0], vertices[2]);
    _mesh.add_edge(vertices[0], vertices[3]);
    _mesh.add_edge(vertices[1], vertices[2]);
    _mesh.add_edge(vertices[2], vertices[3]);
    _mesh.add_edge(vertices[3], vertices[1]);
    _mesh.add_edge(vertices[1], vertices[4]);
    _mesh.add_edge(vertices[2], vertices[4]);
    _mesh.add_edge(vertices[3], vertices[4]);
    _mesh.add_edge(vertices[0], vertices[4]);

}

void TetrahedralMeshBase::generate_tri_tet_with_faces(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle  vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<std::vector<VertexHandle>> faces;
    faces.push_back({vertices[0], vertices[1], vertices[2]});
    faces.push_back({vertices[0], vertices[1], vertices[3]});
    faces.push_back({vertices[0], vertices[2], vertices[3]});
    faces.push_back({vertices[1], vertices[2], vertices[4]});
    faces.push_back({vertices[1], vertices[3], vertices[4]});
    faces.push_back({vertices[2], vertices[3], vertices[4]});
    faces.push_back({vertices[0], vertices[1], vertices[4]});
    faces.push_back({vertices[0], vertices[2], vertices[4]});
    faces.push_back({vertices[0], vertices[3], vertices[4]});

    for (auto face : faces) {
        _mesh.add_face(face);
    }
}

//TODO: this is no tetmesh anymore, as there are polyhedra with 5 faces, 2 triangles and 3 quadrangles
void TetrahedralMeshBase::generate_nested_tets(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[8];
    for (int i = 0; i < 8; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<std::vector<VertexHandle>> faces;
    // inner tet
    faces.push_back({vertices[0], vertices[1], vertices[2]});
    faces.push_back({vertices[0], vertices[1], vertices[3]});
    faces.push_back({vertices[0], vertices[2], vertices[3]});
    faces.push_back({vertices[1], vertices[2], vertices[3]});
    // outer tet
    faces.push_back({vertices[4], vertices[5], vertices[6]});
    faces.push_back({vertices[4], vertices[5], vertices[7]});
    faces.push_back({vertices[4], vertices[6], vertices[7]});
    faces.push_back({vertices[5], vertices[6], vertices[7]});
    // faces between tets
    faces.push_back({vertices[0], vertices[1], vertices[4], vertices[5]});
    faces.push_back({vertices[0], vertices[2], vertices[4], vertices[6]});
    faces.push_back({vertices[1], vertices[2], vertices[5], vertices[6]});
    faces.push_back({vertices[0], vertices[3], vertices[4], vertices[7]});
    faces.push_back({vertices[1], vertices[3], vertices[5], vertices[7]});
    faces.push_back({vertices[2], vertices[3], vertices[6], vertices[7]});

    for (auto face : faces) {
        _mesh.add_face(face);
    }
}

void TetrahedralMeshBase::generate_tet_3F(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 4; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<std::vector<VertexHandle>> faces;
    faces.push_back({vertices[0], vertices[1], vertices[2]});
    faces.push_back({vertices[0], vertices[1], vertices[3]});
    faces.push_back({vertices[0], vertices[2], vertices[3]});

    for (auto face : faces) {
        _mesh.add_face(face);
    }

}

void TetrahedralMeshBase::generate_tets_two_connected_components(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[8];
    for (int i = 0; i < 8; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    std::vector<std::vector<VertexHandle>> cells;
    cells = {{vertices[0], vertices[1], vertices[2], vertices[3]},
             {vertices[4], vertices[5], vertices[6], vertices[7]}};

    for (auto cell : cells) {
        _mesh.add_cell(cell);
    }
}

void TetrahedralMeshBase::generate_non_manifold_tet_3T1F(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[6];
    for (int i = 0; i < 6; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    _mesh.add_edge(vertices[0], vertices[1]);
    _mesh.add_edge(vertices[0], vertices[2]);
    _mesh.add_edge(vertices[0], vertices[3]);
    _mesh.add_edge(vertices[1], vertices[2]);
    _mesh.add_edge(vertices[2], vertices[3]);
    _mesh.add_edge(vertices[3], vertices[1]);
    _mesh.add_edge(vertices[1], vertices[4]);
    _mesh.add_edge(vertices[2], vertices[4]);
    _mesh.add_edge(vertices[3], vertices[4]);
    _mesh.add_edge(vertices[1], vertices[5]);
    _mesh.add_edge(vertices[2], vertices[5]);
    _mesh.add_edge(vertices[3], vertices[5]);
}