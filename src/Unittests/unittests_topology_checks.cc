#include "unittests_common.hh"
#include "OpenVolumeMesh/FileManager/FileManager.hh"

#include <OpenVolumeMesh/Core/TopologyChecks.hh>

using namespace OpenVolumeMesh;

TEST_F(TetrahedralMeshBase, cell_exists) {

    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh); // one simple tet
    auto cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    auto cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    // all four vertices should define the tet
    EXPECT_EQ(cell, *mesh.cells_begin());

    // test, if cell is not found, if one vertex is missing
    generateTetrahedralMesh(mesh);
    cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    cell_vertices.pop_back();
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    EXPECT_EQ(cell, CellHandle(-1));

    // test, if cell is not found for four wrong vertices
    cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    Vec3d p5(0.0, 0.0, 1.0);
    VertexHandle v5 = mesh.add_vertex(p5);
    auto cell_vertices_new = {cell_vertices[0], cell_vertices[1], cell_vertices[2], v5 };
    auto cell_new = mesh.add_cell(cell_vertices_new);
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices_new);
    EXPECT_EQ(cell, cell_new); // just to make sure, the new cell is setup correctly
    auto cell_vertices_inexistant = {cell_vertices[0], cell_vertices[1], cell_vertices[3], v5 };
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices_inexistant);
    EXPECT_EQ(cell, CellHandle(-1));
}

TEST_F(TetrahedralMeshBase, face_contains_vertex) {

    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    auto face = *mesh.faces_begin();
    auto vertex_in_face = *mesh.fv_iter(face);
    EXPECT_TRUE(OpenVolumeMesh::face_contains_vertex(mesh, vertex_in_face, face));
    VH vertex_outside_face = VH(-1);
    for (auto vertex : mesh.vertices()) {
        //TODO: this can probably be done a lot nicer
        bool found_vertex = false;
        for (auto face_vertex : mesh.face_vertices(face)) {
            if (vertex == face_vertex) {
                found_vertex = true;
                continue;
            }
        }
        if (!found_vertex) {
            vertex_outside_face = vertex;
            break;
        }
    }
    EXPECT_FALSE(OpenVolumeMesh::face_contains_vertex(mesh, vertex_outside_face, face));
}

TEST_F(TetrahedralMeshBase, cell_contains_vertex) {

    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    auto cell = *mesh.cells_begin();
    auto vertex_in_cell = *mesh.cv_iter(cell);
    EXPECT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex_in_cell, cell));
    auto cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    Vec3d p5(0.0, 0.0, 1.0);
    VertexHandle v5 = mesh.add_vertex(p5);
    auto cell_vertices_new = {cell_vertices[0], cell_vertices[2], cell_vertices[1], v5 };
    auto new_cell = mesh.add_cell(cell_vertices_new);
    EXPECT_FALSE(OpenVolumeMesh::cell_contains_vertex(mesh, cell_vertices[3], new_cell));
    EXPECT_FALSE(OpenVolumeMesh::cell_contains_vertex(mesh, v5, cell));
}

/*
TEST_F(TetrahedralMeshBase, find_non_cell_tets) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);
    std::set<std::pair<std::set<VertexHandle>, bool>> non_cell_tets;

    auto found_non_cell_tets = OpenVolumeMesh::findNonCellTets(mesh):
    EXPECT_TRUE(found_non_cell_tets.empty());
}
 */

TEST_F(TetrahedralMeshBase, singleConnectedComponent) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::singleConnectedComponent(mesh));

    Vec3d p5(0.0, 0.0, 1.0);
    VertexHandle v5 = mesh.add_vertex(p5);
    EXPECT_FALSE(OpenVolumeMesh::singleConnectedComponent(mesh));
}

TEST_F(TetrahedralMeshBase, containsVoid) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    EXPECT_FALSE(OpenVolumeMesh::containsVoid(mesh));

    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.readFile("Cuboid.ovm", mesh);

    CellHandle cell_to_remove(-1);
    for (auto cell : mesh.cells()){
        if (mesh.is_boundary(cell)) continue;
        bool has_boundary_vertex = false;
        for (auto cv_it = mesh.cv_iter(cell); cv_it.is_valid(); ++cv_it) {
            if (mesh.is_boundary(*cv_it)) {
                has_boundary_vertex = true;
                break;
            }
        }
        if (!has_boundary_vertex) {
            cell_to_remove = cell;
            break;
        }
    }

    ASSERT_TRUE(cell_to_remove.is_valid());

    mesh.delete_cell(cell_to_remove);
    EXPECT_TRUE(OpenVolumeMesh::containsVoid(mesh));
}

TEST_F(TetrahedralMeshBase, manifoldVertex) {
    // test some manifold meshes
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    auto vertex = *mesh.vertices_begin();

    EXPECT_TRUE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateTetrahedralMesh_2(mesh);
    //Assume vertex 0 is in the common edge between the two cells
    vertex = *mesh.vertices_begin();
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));

    EXPECT_TRUE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    // test some non-manifold meshes
    generateNonManifoldTet_2T1V(mesh);
    // Assume vertex 0 is the common vertex between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *++mesh.cells_begin()));

    vertex = *mesh.vertices_begin();

    for (CellHandle c : mesh.vertex_cells(vertex))
        std::cout << c << std::endl;

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_2T1E(mesh);
    vertex = *mesh.vertices_begin();
    auto vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *++mesh.cells_begin()));

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_6T1E(mesh);
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two fans
    for (auto c_it = mesh.cells_begin(); c_it.is_valid(); ++c_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *c_it));
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *c_it));
    }

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_1T1F1E(mesh);
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the tet and the triangle
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_EQ(mesh.valence(vertex), 4);
    ASSERT_EQ(mesh.valence(vertex2), 4);

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_2F1E(mesh);
    vertex = *mesh.vertices_begin();

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_3T1V3E(mesh);
    vertex = *mesh.vertices_begin();
    // Assume vertex 0 is shared among all three tets
    int count = 0;
    for (auto vc_it = mesh.vc_iter(vertex); vc_it.is_valid(); ++vc_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *vc_it));
        ++count;
    }
    ASSERT_EQ(count, 3);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1V(mesh);
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1E(mesh);
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1F(mesh);
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1T1F1V(mesh);
    // Assume vertex 0 is shared by the tet and the face
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 5);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_2F1V(mesh);
    // Assume vertex 0 is shared by the two faces
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 4);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1F1E1V(mesh);
    // Assume vertex 0 is share by the face and the edge
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 3);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
}

TEST_F(TetrahedralMeshBase, noDoubleEdges) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::noDoubleEdges(mesh));

    auto heh = *mesh.halfedges_begin();
    auto from_vertex = mesh.from_vertex_handle(heh);
    auto to_vertex = mesh.to_vertex_handle(heh);
    mesh.add_edge(from_vertex, to_vertex, true);

    EXPECT_FALSE(OpenVolumeMesh::noDoubleEdges(mesh));
}