#include "unittests_common.hh"
#include "OpenVolumeMesh/FileManager/FileManager.hh"

#include <OpenVolumeMesh/Core/TopologyChecks.hh>

using namespace OpenVolumeMesh;

TEST_F(TetrahedralMeshBase, cell_exists) {

    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh); // one simple tet
    auto cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    auto cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    // all four vertices should define the tet
    EXPECT_EQ(cell, *mesh.cells_begin());

    // test, if cell is not found, if one vertex is missing
    generate_tetrahedral_mesh(mesh);
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
    generate_tetrahedral_mesh(mesh);
    std::vector<VertexHandle> empty_vertices;
    cell = OpenVolumeMesh::cell_exists(mesh, empty_vertices);
    EXPECT_EQ(cell, CellHandle(-1)); // valid mesh but invalid vertex list
    cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());
    mesh.clear();
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    EXPECT_EQ(cell, CellHandle(-1)); // empty mesh but valid vertex list
    cell = OpenVolumeMesh::cell_exists(mesh, empty_vertices);
    EXPECT_EQ(cell, CellHandle(-1)); // empty mesh and empty vertex list
}

TEST_F(TetrahedralMeshBase, face_contains_vertex) {

    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

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
    generate_tetrahedral_mesh(mesh);

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

TEST_F(TetrahedralMeshBase, single_connected_component) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::single_connected_component(mesh));

    Vec3d p5(0.0, 0.0, 1.0);
    VertexHandle v5 = mesh.add_vertex(p5);
    EXPECT_FALSE(OpenVolumeMesh::single_connected_component(mesh));

    // test invalid input
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::single_connected_component(mesh));
}

TEST_F(TetrahedralMeshBase, contains_void) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_FALSE(OpenVolumeMesh::contains_void(mesh));

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
    EXPECT_TRUE(OpenVolumeMesh::contains_void(mesh));

    // Test invalid inputs
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::contains_void(mesh));
}

TEST_F(TetrahedralMeshBase, manifold_vertex) {
    // ------------------------------ //
    // test some manifold meshes
    // ------------------------------ //
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh); // One simple tet

    auto vertex = *mesh.vertices_begin();

    EXPECT_TRUE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_tetrahedral_mesh_2(mesh); // Two tets sharing exactly one face
    //Assume vertex 0 is in the common edge between the two cells
    vertex = *mesh.vertices_begin();
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));

    EXPECT_TRUE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    // ------------------------------ //
    // test some non-manifold meshes
    // ------------------------------ //
    generate_non_manifold_tet_2T1V(mesh); // Two tets sharing exactly one vertex
    // Assume vertex 0 is the common vertex between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, *mesh.vertices_begin(), *++mesh.cells_begin()));

    vertex = *mesh.vertices_begin();

    for (CellHandle c : mesh.vertex_cells(vertex))
        std::cout << c << std::endl;

    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_2T1E(mesh); // Two tets sharing exactly one edge
    vertex = *mesh.vertices_begin();
    auto vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two cells
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *++mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *++mesh.cells_begin()));

    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex2));

    generate_non_manifold_tet_6T1E(mesh); // Two fans of three tets each sharing exactly one edge
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the two fans
    for (auto c_it = mesh.cells_begin(); c_it.is_valid(); ++c_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *c_it));
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *c_it));
    }

    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex2));

    generate_non_manifold_tet_1T1F1E(mesh); // One tet and one face sharing exactly one edge
    vertex = *mesh.vertices_begin();
    vertex2 = *++mesh.vertices_begin();
    // Assume vertices 0 and 1 are in the common edge between the tet and the triangle
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *mesh.cells_begin()));
    ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex2, *mesh.cells_begin()));
    ASSERT_EQ(mesh.valence(vertex), 4);
    ASSERT_EQ(mesh.valence(vertex2), 4);

    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex2));

    generate_non_manifold_tet_2F1E(mesh); // Two faces sharing exactly one edge
    vertex = *mesh.vertices_begin();

    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_3T1V3E(mesh); // Three tets sharing exactly one vertex and pairwise sharing exactly one edge
    vertex = *mesh.vertices_begin();
    // Assume vertex 0 is shared among all three tets
    int count = 0;
    for (auto vc_it = mesh.vc_iter(vertex); vc_it.is_valid(); ++vc_it) {
        ASSERT_TRUE(OpenVolumeMesh::cell_contains_vertex(mesh, vertex, *vc_it));
        ++count;
    }
    ASSERT_EQ(count, 3);
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_1V(mesh); // One isolated vertex
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_1E(mesh); // One isolated edge
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_1F(mesh); // One isolated face
    vertex = *mesh.vertices_begin();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_1T1F1V(mesh); // One tet and one face sharing exactly one vertex
    // Assume vertex 0 is shared by the tet and the face
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 5);
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_2F1V(mesh); // Two faces sharing exactly one vertex
    // Assume vertex 0 is shared by the two faces
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 4);
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    generate_non_manifold_tet_1F1E1V(mesh); //  One tet and one edge sharing exactly one vertex
    // Assume vertex 0 is share by the face and the edge
    vertex = *mesh.vertices_begin();
    ASSERT_EQ(mesh.valence(vertex), 3);
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex));

    // ------------------------------ //
    // Test some invalid inputs
    // ------------------------------ //
    vertex.reset();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex)); // vertex not valid
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex)); // mesh and vertex not valid
    generate_non_manifold_tet_1V(mesh); // One isolated vertex
    vertex = *mesh.vertices_begin();
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex)); // mesh not valid
    generate_non_manifold_tet_1V(mesh);
    EXPECT_FALSE(OpenVolumeMesh::manifold_vertex(mesh, vertex)); // mesh and vertex valid, but vertex not in mesh
}

TEST_F(TetrahedralMeshBase, no_double_edges) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::no_double_edges(mesh));

    auto heh = *mesh.halfedges_begin();
    auto from_vertex = mesh.from_vertex_handle(heh);
    auto to_vertex = mesh.to_vertex_handle(heh);
    mesh.add_edge(from_vertex, to_vertex, true);

    EXPECT_FALSE(OpenVolumeMesh::no_double_edges(mesh));

    // Test invalid input
    mesh.clear();
    EXPECT_TRUE(OpenVolumeMesh::no_double_edges(mesh));
}

TEST_F(TetrahedralMeshBase, find_non_cell_tets) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::find_non_cell_tets(mesh, true).empty());


    generate_non_manifold_tet_3T1V3E(mesh);

    // find the "missing tet" in the mesh manually
    std::set<VertexHandle> tet_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        // there are 7 vertices. 3 with valence 3 outside the missing tet, 3 with valence 5, defining one face of the
        // missing tet and one central vertex with valence 6, being the last vertex of the missing tet.
        if (mesh.valence(*v_it) == 5 || mesh.valence(*v_it) == 6) {
            tet_vertices.insert(*v_it);
        }
    }
    ASSERT_EQ(tet_vertices.size(), 4);
    std::vector<VertexHandle> face_vertices = {*++mesh.vertices_begin(), *++(++(++mesh.vertices_begin())), *++(++mesh.vertices_begin())};
//    print_mesh_topology(mesh);
    mesh.add_face(face_vertices);
//    print_mesh_topology(mesh);
    std::set<std::set<VertexHandle>> res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 1);
    std::set<VertexHandle> found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);

    generate_tet_without_cell(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 1);
    tet_vertices.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        tet_vertices.insert(*v_it);
    }
    found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);

    generate_tri_tet_with_faces(mesh);
    print_mesh_topology(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 3);
    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::vector<std::vector<VertexHandle>> expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[2], all_vertices[3], all_vertices[4]}
    };
    std::set<std::set<VertexHandle>> expected_tets;
    std::set<VertexHandle> exptected_tet;
    for (auto expected_tet_vertices : expected_tets_vertices) {
        exptected_tet.clear();
        for (auto expected_vertex : expected_tet_vertices) {
            exptected_tet.insert(expected_vertex);
        }
        expected_tets.insert(exptected_tet);
    }
    EXPECT_EQ(expected_tets, res);

    generate_nested_tets(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 2);
//    std::vector<VertexHandle> all_vertices;
    all_vertices.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
//    std::vector<std::vector<VertexHandle>> expected_tets_vertices = {
    expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[3]},
            {all_vertices[4], all_vertices[5], all_vertices[6], all_vertices[7]}
    };
//    std::set<std::set<VertexHandle>> expected_tets;
    expected_tets.clear();
//    std::set<VertexHandle> exptected_tet;
    for (auto expectedTetVertices : expected_tets_vertices) {
        exptected_tet.clear();
        for (auto expected_vertex : expectedTetVertices) {
            exptected_tet.insert(expected_vertex);
        }
        expected_tets.insert(exptected_tet);
    }
    EXPECT_EQ(expected_tets, res);

    // Don't detect a tet with some faces missing
    generate_tet_3F(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 0);

    // test some invalid inputs
    mesh.clear();
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(res.empty());
}

TEST_F(TetrahedralMeshBase, find_non_face_triangles_around_vertex) {
    TetrahedralMesh &mesh = this->mesh_;

    // find the non-face triangle
    generate_tri_without_face(mesh);
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
    generate_non_manifold_tet_1F(mesh);
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
    generate_tri_without_face(mesh);
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
    generate_non_manifold_tet_1F(mesh);
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);

    EXPECT_EQ(found_tris.size(), 0);

    // find non-face triangles in a mesh consisting of two tets without cells or faces
    generate_tet_without_cells_and_faces(mesh);
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);
    EXPECT_EQ(found_tris.size(), 7);

    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::set<std::set<VertexHandle>> expected_tris;
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[3]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[0]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[4]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[1]);
    expected_tri.insert(all_vertices[3]);
    expected_tri.insert(all_vertices[4]);
    expected_tris.insert(expected_tri);
    expected_tri.clear();
    expected_tri.insert(all_vertices[2]);
    expected_tri.insert(all_vertices[3]);
    expected_tri.insert(all_vertices[4]);
    expected_tris.insert(expected_tri);
    EXPECT_EQ(expected_tris, found_tris);

    // test some invalid inputs
    mesh.clear();
    found_tris = OpenVolumeMesh::find_non_face_triangles(mesh);
    EXPECT_TRUE(found_tris.empty());
}

TEST_F(TetrahedralMeshBase, findNonCellTets_nonfaceTris) {
    TetrahedralMesh &mesh = this->mesh_;

    generate_tet_without_cells_and_faces(mesh);

    // find expected tets
    std::vector<VertexHandle> all_vertices;
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    std::set<std::set<VertexHandle>> expected_tets;
    std::set<VertexHandle> expected_tet;
    expected_tet.insert(all_vertices[0]);
    expected_tet.insert(all_vertices[1]);
    expected_tet.insert(all_vertices[2]);
    expected_tet.insert(all_vertices[3]);
    expected_tets.insert(expected_tet);
    expected_tet.clear();
    expected_tet.insert(all_vertices[1]);
    expected_tet.insert(all_vertices[2]);
    expected_tet.insert(all_vertices[3]);
    expected_tet.insert(all_vertices[4]);
    expected_tets.insert(expected_tet);

    auto found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(found_tets.size(), 2);
    EXPECT_EQ(found_tets, expected_tets);

    // test mixed case
    std::vector<VertexHandle> vertices_cell0(expected_tets.begin()->size());
    std::copy(expected_tets.begin()->begin(), expected_tets.begin()->end(), vertices_cell0.begin());
    mesh.add_cell(vertices_cell0);
    expected_tets.erase(expected_tets.begin());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(found_tets.size(), 1);
    EXPECT_EQ(found_tets, expected_tets);

    generate_non_manifold_tet_3T1F(mesh);
    all_vertices.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }
    expected_tets.clear();
    expected_tet.clear();
    expected_tet.insert(all_vertices[0]);
    expected_tet.insert(all_vertices[1]);
    expected_tet.insert(all_vertices[2]);
    expected_tet.insert(all_vertices[3]);
    expected_tets.insert(expected_tet);
    expected_tet.clear();
    expected_tet.insert(all_vertices[1]);
    expected_tet.insert(all_vertices[2]);
    expected_tet.insert(all_vertices[3]);
    expected_tet.insert(all_vertices[4]);
    expected_tets.insert(expected_tet);
    expected_tet.clear();
    expected_tet.insert(all_vertices[1]);
    expected_tet.insert(all_vertices[2]);
    expected_tet.insert(all_vertices[3]);
    expected_tet.insert(all_vertices[5]);
    expected_tets.insert(expected_tet);

    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(found_tets.size(), 3);
    EXPECT_EQ(found_tets, expected_tets);

    // check that no existing tet is found
    generate_tetrahedral_mesh(mesh);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_TRUE(found_tets.empty());

    // test some invalid inputs
    mesh.clear();
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_TRUE(found_tets.empty());
}