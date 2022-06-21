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

void TetrahedralMeshBase::generateTetrahedralMesh(TetrahedralMesh& _mesh) {

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

void TetrahedralMeshBase::generateTetrahedralMesh_2(TetrahedralMesh &_mesh) {
    _mesh.clear();

    VertexHandle vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = _mesh.add_vertex();
    }
    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3], true);
    _mesh.add_cell(vertices[0],  vertices[2], vertices[1], vertices[4]);
}

void TetrahedralMeshBase::generateNonManifoldTet_2T1V(TetrahedralMesh& mesh) {

    mesh.clear();
/*
    Vec3d positions[7] = {
            Vec3d(0.0, 0.0, 0.0),
            Vec3d(1.0, 0.0, 0.0),
            Vec3d(0.0, 1.0, 0.0),
            Vec3d(0.0, 0.0, 1.0),
            Vec3d(-1.0, 0.0, 0.0),
            Vec3d(0.0, -1.0, 0.0),
            Vec3d(0.0, 0.0, -1.0)
    };
    */
    VertexHandle vertices[7];
    for (int i = 0; i < 7; ++i) {
//        vertices[i] = mesh.add_vertex(positions[i]);
        vertices[i] = mesh.add_vertex();
    }
    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    mesh.add_cell(vertices[0], vertices[4], vertices[5], vertices[6]);
}

void TetrahedralMeshBase::generateNonManifoldTet_2T1E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[6];
    for (int i = 0; i < 6; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    mesh.add_cell(vertices[0], vertices[1], vertices[4], vertices[5]);
}

void TetrahedralMeshBase::generateNonManifoldTet_6T1E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[8];
    for (int i = 0; i < 8; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    mesh.add_cell(vertices[1], vertices[0], vertices[2], vertices[4]);
    mesh.add_cell(vertices[0], vertices[1], vertices[5], vertices[6]);
    mesh.add_cell(vertices[1], vertices[0], vertices[5], vertices[7]);
}

void TetrahedralMeshBase::generateNonManifoldTet_1T1F1E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3], true);
    std::vector<VertexHandle> face_vertices = {vertices[0], vertices[1], vertices[4]};
    mesh.add_face(face_vertices);
}

void TetrahedralMeshBase::generateNonManifoldTet_2F1E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    std::vector<VertexHandle> faceVertices = {vertices[0], vertices[1], vertices[2]};
    mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generateNonManifoldTet_3T1V3E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle  vertices[7];
    for (int i = 0; i < 7; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[4]);
    mesh.add_cell(vertices[0], vertices[1], vertices[3], vertices[5]);
    mesh.add_cell(vertices[0], vertices[2], vertices[3], vertices[6]);
}

void TetrahedralMeshBase::generateNonManifoldTet_1V(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle v = mesh.add_vertex();
}

void TetrahedralMeshBase::generateNonManifoldTet_1E(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle  vertices[2];
    for (int i = 0; i < 2; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_edge(vertices[0], vertices[1]);
}

void TetrahedralMeshBase::generateNonManifoldTet_1F(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[3];
    for (int i = 0; i < 3; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    std::vector<VertexHandle> faceVertices = {vertices[0], vertices[1], vertices[2]};
    mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generateNonManifoldTet_1T1F1V(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle  vertices[6];
    for (int i = 0; i < 6; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    std::vector<VertexHandle> faceVertices = { vertices[0], vertices[4], vertices[5]};
    mesh.add_face(faceVertices);
}

void TetrahedralMeshBase::generateNonManifoldTet_2F1V(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle  vertices[5];
    for (int i = 0; i < 5; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    std::vector<VertexHandle> verticesFace1 = {vertices[0], vertices[1], vertices[2]};
    std::vector<VertexHandle> verticesFace2 = {vertices[0], vertices[3], vertices[4]};
    mesh.add_face(verticesFace1);
    mesh.add_face(verticesFace2);
}

void TetrahedralMeshBase::generateNonManifoldTet_1F1E1V(TetrahedralMesh &mesh) {
    mesh.clear();

    VertexHandle vertices[4];
    for (int i = 0; i < 4; ++i) {
        vertices[i] = mesh.add_vertex();
    }

    std::vector<VertexHandle> faceVertices = {vertices[0], vertices[1], vertices[2]};
    mesh.add_face(faceVertices);
    mesh.add_edge(vertices[0], vertices[3]);
}