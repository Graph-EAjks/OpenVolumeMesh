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

    // test, if it returns an invalid cell for an invalid input
    generateTetrahedralMesh(mesh);
    std::vector<VertexHandle> emptyVertices;
    cell = OpenVolumeMesh::cell_exists(mesh, emptyVertices);
    EXPECT_EQ(cell, CellHandle(-1)); // valid mesh but invalid vertex list
    cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    mesh.clear();
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    EXPECT_EQ(cell, CellHandle(-1)); // empty mesh but valid vertex list
    cell = OpenVolumeMesh::cell_exists(mesh, emptyVertices);
    EXPECT_EQ(cell, CellHandle(-1)); // empty mesh and empty vertex list
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

    // Test some invalid inputs
    vertex_outside_face = VertexHandle(-1);
    EXPECT_FALSE(OpenVolumeMesh::face_contains_vertex(mesh, vertex_outside_face, face)); // invalid vertex
    face = FaceHandle(-1);
    EXPECT_FALSE(OpenVolumeMesh::face_contains_vertex(mesh, vertex_in_face, face)); // invalid face
    face = *mesh.faces_begin();
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::face_contains_vertex(mesh, vertex_in_face, face)); // invalid mesh
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

    // Test some invalid inputs
    v5 = VertexHandle(-1);
    EXPECT_FALSE(OpenVolumeMesh::cell_contains_vertex(mesh, v5, cell)); // invalid vertex
    cell = CellHandle(-1);
    EXPECT_FALSE(OpenVolumeMesh::cell_contains_vertex(mesh, cell_vertices[0], cell)); // invalid cell
    cell = *mesh.cells_begin();
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::cell_contains_vertex(mesh, cell_vertices[0], cell)); // invalid mesh
}

TEST_F(TetrahedralMeshBase, singleConnectedComponent) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::singleConnectedComponent(mesh));

    Vec3d p5(0.0, 0.0, 1.0);
    VertexHandle v5 = mesh.add_vertex(p5);
    EXPECT_FALSE(OpenVolumeMesh::singleConnectedComponent(mesh));

    // test invalid input
    mesh.clear();
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

    // Test invalid inputs
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::containsVoid(mesh));
}

TEST_F(TetrahedralMeshBase, manifoldVertex) {
    // ------------------------------ //
    // test some manifold meshes
    // ------------------------------ //
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh); // One simple tet

    auto vertex = *mesh.vertices_begin();

    EXPECT_TRUE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateTetrahedralMesh_2(mesh); // Two tets sharing exactly one face
    //Assume vertex 0 is in the common edge between the two cells
    vertex = *mesh.vertices_begin();
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));

    EXPECT_TRUE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    // ------------------------------ //
    // test some non-manifold meshes
    // ------------------------------ //
    generateNonManifoldTet_2T1V(mesh); // Two tets sharing exactly one vertex
    // Assume vertex 0 is the common vertex between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *++mesh.cells_begin()));

    vertex = *mesh.vertices_begin();

    for (CellHandle c : mesh.vertex_cells(vertex))
        std::cout << c << std::endl;

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_2T1E(mesh); // Two tets sharing exactly one edge
    vertex = *mesh.vertices_begin();
    auto vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *++mesh.cells_begin()));

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_6T1E(mesh); // Two fans of three tets each sharing exactly one edge
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two fans
    for (auto c_it = mesh.cells_begin(); c_it.is_valid(); ++c_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *c_it));
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *c_it));
    }

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_1T1F1E(mesh); // One tet and one face sharing exactly one edge
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the tet and the triangle
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_EQ(mesh.valence(vertex), 4);
    ASSERT_EQ(mesh.valence(vertex2), 4);

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex2));

    generateNonManifoldTet_2F1E(mesh); // Two faces sharing exactly one edge
    vertex = *mesh.vertices_begin();

    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_3T1V3E(mesh); // Three tets sharing exactly one vertex and pairwise sharing exactly one edge
    vertex = *mesh.vertices_begin();
    // Assume vertex 0 is shared among all three tets
    int count = 0;
    for (auto vc_it = mesh.vc_iter(vertex); vc_it.is_valid(); ++vc_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *vc_it));
        ++count;
    }
    ASSERT_EQ(count, 3);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1V(mesh); // One isolated vertex
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1E(mesh); // One isolated edge
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1F(mesh); // One isolated face
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1T1F1V(mesh); // One tet and one face sharing exactly one vertex
    // Assume vertex 0 is shared by the tet and the face
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 5);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_2F1V(mesh); // Two faces sharing exactly one vertex
    // Assume vertex 0 is shared by the two faces
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 4);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    generateNonManifoldTet_1F1E1V(mesh); //  One tet and one edge sharing exactly one vertex
    // Assume vertex 0 is share by the face and the edge
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 3);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex));

    // ------------------------------ //
    // Test some invalid inputs
    // ------------------------------ //
    vertex.reset();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex)); // vertex not valid
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex)); // mesh and vertex not valid
    generateNonManifoldTet_1V(mesh); // One isolated vertex
    vertex = *mesh.vertices_begin();
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex)); // mesh not valid
    generateNonManifoldTet_1V(mesh);
    EXPECT_FALSE(OpenVolumeMesh::manifoldVertex(mesh, vertex)); // mesh and vertex valid, but vertex not in mesh
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

    // Test invalid input
    mesh.clear();
    EXPECT_TRUE(OpenVolumeMesh::noDoubleEdges(mesh));
}

TEST_F(TetrahedralMeshBase, findNonCellTets) {
    TetrahedralMesh &mesh = this->mesh_;
    generateTetrahedralMesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::findNonCellTets(mesh, true).empty());


    generateNonManifoldTet_3T1V3E(mesh);

    // find the "missing tet" in the mesh manually
    std::set<VertexHandle> tetVertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        // there are 7 vertices. 3 with valence 3 outside the missing tet, 3 with valence 5, defining one face of the
        // missing tet and one central vertex with valence 6, being the last vertex of the missing tet.
        if (mesh.valence(*v_it) == 5 || mesh.valence(*v_it) == 6) {
            tetVertices.insert(*v_it);
        }
    }
    ASSERT_EQ(tetVertices.size(), 4);
    std::vector<VertexHandle> face_vertices = {*++mesh.vertices_begin(), *++(++(++mesh.vertices_begin())), *++(++mesh.vertices_begin())};
//    printMeshTopology(mesh);
    mesh.add_face(face_vertices);
//    printMeshTopology(mesh);
    std::set<std::set<VertexHandle>> res = OpenVolumeMesh::findNonCellTets(mesh, true);
    EXPECT_EQ(res.size(), 1);
    std::set<VertexHandle> foundTet = *res.begin();
    EXPECT_EQ(tetVertices, foundTet);

    generateTetWithoutCell(mesh);
    res = OpenVolumeMesh::findNonCellTets(mesh, true);
    EXPECT_EQ(res.size(), 1);
    tetVertices.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        tetVertices.insert(*v_it);
    }
    foundTet = *res.begin();
    EXPECT_EQ(tetVertices, foundTet);

    generateTriTet_withFaces(mesh);
    printMeshTopology(mesh);
    res = OpenVolumeMesh::findNonCellTets(mesh, true);
    EXPECT_EQ(res.size(), 3);
    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::vector<std::vector<VertexHandle>> expectedTetsVertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[2], all_vertices[3], all_vertices[4]}
    };
    std::set<std::set<VertexHandle>> expectedTets;
    std::set<VertexHandle> exptectedTet;
    for (auto expectedTetVertices : expectedTetsVertices) {
        exptectedTet.clear();
        for (auto expectedVertex : expectedTetVertices) {
            exptectedTet.insert(expectedVertex);
        }
        expectedTets.insert(exptectedTet);
    }
    EXPECT_EQ(expectedTets, res);

    // test some invalid inputs
    mesh.clear();
    res = OpenVolumeMesh::findNonCellTets(mesh, true);
    EXPECT_TRUE(res.empty());
}

TEST_F(TetrahedralMeshBase, find_non_face_triangles_around_vertex) {
    TetrahedralMesh &mesh = this->mesh_;

    // find the non-face triangle
    generateTriWithoutFace(mesh);
    VertexHandle vertex = *mesh.vertices_begin();
    auto found_tris = OpenVolumeMesh::find_non_face_triangles_around_vertex(mesh, vertex);

    EXPECT_EQ(found_tris.size(), 1);
    std::set<VertexHandle> found_tri = *found_tris.begin();

    EXPECT_EQ(found_tri.size(), 3);
    std::set<VertexHandle> expected_tri;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        expected_tri.insert(*v_it);
    }
    std::set<VertexHandle> found_tri_set;
    for (auto v : found_tri) {
        found_tri_set.insert(v);
    }

    EXPECT_EQ(expected_tri, found_tri_set);

    // Don't detect any existing faces
    generateNonManifoldTet_1F(mesh);
    vertex = *mesh.vertices_begin();
    found_tris = OpenVolumeMesh::find_non_face_triangles_around_vertex(mesh, vertex);

    EXPECT_EQ(found_tris.size(), 0);

    // test some invalid inputs
    vertex = VertexHandle(-1);
    found_tris = OpenVolumeMesh::find_non_face_triangles_around_vertex(mesh, vertex);
    EXPECT_TRUE(found_tris.empty());

    vertex = *mesh.vertices_begin();
    mesh.clear();
    found_tris = OpenVolumeMesh::find_non_face_triangles_around_vertex(mesh, vertex);
    EXPECT_TRUE(found_tris.empty());

    vertex = VertexHandle(-1);
    found_tris = OpenVolumeMesh::find_non_face_triangles_around_vertex(mesh, vertex);
    EXPECT_TRUE(found_tris.empty());
}

TEST_F(TetrahedralMeshBase, find_non_face_triangles) {
    TetrahedralMesh &mesh = this->mesh_;

    // same tests as for the find_non_face_triangles_around_vertex test
    // find the non-face triangle
    generateTriWithoutFace(mesh);
    auto found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);

    EXPECT_EQ(found_tris.size(), 1);
    std::set<VertexHandle> found_tri = *found_tris.begin();

    EXPECT_EQ(found_tri.size(), 3);
    std::set<VertexHandle> expected_tri;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        expected_tri.insert(*v_it);
    }
    std::set<VertexHandle> found_tri_set;
    for (auto v : found_tri) {
        found_tri_set.insert(v);
    }

    EXPECT_EQ(expected_tri, found_tri_set);

    // Don't detect any existing faces
    generateNonManifoldTet_1F(mesh);
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);

    EXPECT_EQ(found_tris.size(), 0);

    // find non-face triangles in a mesh consisting of two tets without cells or faces
    generateTet_withoutCellsAndFaces(mesh);
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);
    EXPECT_EQ(found_tris.size(), 7);

    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::set<std::set<VertexHandle>> expectedTris;
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[3]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[4]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[3]);
    expected_tri.insert(all_vertices[4]);
    expectedTris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expected_tri.insert(all_vertices[4]);
    expectedTris.insert(expected_tri);
    EXPECT_EQ(expectedTris, found_tris);

    // test some invalid inputs
    mesh.clear();
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);
    EXPECT_TRUE(found_tris.empty());
}

TEST_F(TetrahedralMeshBase, findNonCellTets_nonfaceTris) {
    TetrahedralMesh &mesh = this->mesh_;

    generateTet_withoutCellsAndFaces(mesh);

    // find expected tets
    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::set<std::set<VertexHandle>> expectedTets;
    std::set<VertexHandle> expectedTet;
    expectedTet.insert(all_vertices[0]);
    expectedTet.insert(all_vertices[1]);
    expectedTet.insert(all_vertices[2]);
    expectedTet.insert(all_vertices[3]);
    expectedTets.insert(expectedTet);
    expectedTet.clear();
    expectedTet.insert(all_vertices[1]);
    expectedTet.insert(all_vertices[2]);
    expectedTet.insert(all_vertices[3]);
    expectedTet.insert(all_vertices[4]);
    expectedTets.insert(expectedTet);

    auto foundTets = OpenVolumeMesh::findNonCellTets(mesh, false);
    EXPECT_EQ(foundTets.size(), 2);
    EXPECT_EQ(foundTets, expectedTets);

    // test mixed case
    std::vector<VertexHandle> vertices_cell0(expectedTets.begin()->size());
    std::copy(expectedTets.begin()->begin(), expectedTets.begin()->end(), vertices_cell0.begin());
    mesh.add_cell(vertices_cell0);
    expectedTets.erase(expectedTets.begin());
    foundTets = OpenVolumeMesh::findNonCellTets(mesh, false);
    EXPECT_EQ(foundTets.size(), 1);
    EXPECT_EQ(foundTets, expectedTets);

    // check that no existing tet is found
    generateTetrahedralMesh(mesh);
    foundTets = OpenVolumeMesh::findNonCellTets(mesh, false);
    EXPECT_TRUE(foundTets.empty());

    // test some invalid inputs
    mesh.clear();
    foundTets = OpenVolumeMesh::findNonCellTets(mesh, false);
    EXPECT_TRUE(foundTets.empty());
}