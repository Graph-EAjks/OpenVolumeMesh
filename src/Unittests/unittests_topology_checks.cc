#include "unittests_common.hh"
#include "OpenVolumeMesh/FileManager/FileManager.hh"

#include <OpenVolumeMesh/Core/TopologyChecks.hh>
#include <OpenVolumeMesh/Core/detail/TopologicalLinkT_impl.hh>

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

    generate_tets_two_connected_components(mesh);
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

TEST_F(TetrahedralMeshBase, contains_double_edges) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_FALSE(OpenVolumeMesh::contains_double_edges(mesh).has_value());

    auto heh = *mesh.halfedges_begin();
    auto from_vertex = mesh.from_vertex_handle(heh);
    auto to_vertex = mesh.to_vertex_handle(heh);
    mesh.add_edge(from_vertex, to_vertex, true);

    auto ret = OpenVolumeMesh::contains_double_edges(mesh);
    EXPECT_TRUE(ret.has_value());
    auto heh_0 = ret->first;
    auto heh_1 = ret->second;
    EXPECT_TRUE((mesh.to_vertex_handle(heh_0) == mesh.to_vertex_handle(heh_1)) ||
                        (mesh.to_vertex_handle(heh_0) == mesh.from_vertex_handle(heh_1)));
    EXPECT_TRUE((mesh.from_vertex_handle(heh_0) == mesh.to_vertex_handle(heh_1)) ||
                (mesh.from_vertex_handle(heh_0) == mesh.from_vertex_handle(heh_1)));

    // Test invalid input
    mesh.clear();
    EXPECT_FALSE(OpenVolumeMesh::contains_double_edges(mesh).has_value());
}

TEST_F(TetrahedralMeshBase, find_multi_edges) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::find_multi_edges(mesh).empty());

    auto heh = *mesh.halfedges_begin();
    auto from_vertex = mesh.from_vertex_handle(heh);
    auto to_vertex = mesh.to_vertex_handle(heh);
    auto eh_2 = mesh.add_edge(from_vertex, to_vertex, true);
    auto heh_2 = mesh.halfedge_handle(eh_2, 0);
    if (mesh.to_vertex_handle(heh_2) != to_vertex) {
        heh_2 = mesh.opposite_halfedge_handle(heh_2);
        ASSERT_EQ(mesh.to_vertex_handle(heh_2), to_vertex);
    }

    EXPECT_EQ(OpenVolumeMesh::find_multi_edges(mesh).size(), 2);
    std::set<std::set<HEH>> expectedMultiEdges;
    std::set<HEH> expectedMultiEdge = {heh, heh_2};
    expectedMultiEdges.insert(expectedMultiEdge);
    expectedMultiEdge = {mesh.opposite_halfedge_handle(heh), mesh.opposite_halfedge_handle(heh_2)};
    expectedMultiEdges.insert(expectedMultiEdge);
    EXPECT_EQ(expectedMultiEdges, OpenVolumeMesh::find_multi_edges(mesh));

    // Test invalid input
    mesh.clear();
    EXPECT_TRUE(OpenVolumeMesh::find_multi_edges(mesh).empty());
}

TEST_F(TetrahedralMeshBase, find_non_cell_tets) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::find_non_cell_tets(mesh, true).empty());
    EXPECT_TRUE(OpenVolumeMesh::find_non_cell_tets(mesh, false).empty());

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
    std::set<std::set<VertexHandle>> res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(res.empty());
    res = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(res.size(), 1);
    std::set<VertexHandle> found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);
    std::vector<VertexHandle> face_vertices = {*++mesh.vertices_begin(), *++(++(++mesh.vertices_begin())), *++(++mesh.vertices_begin())};
//    print_mesh_topology(mesh);
    mesh.add_face(face_vertices);
//    print_mesh_topology(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 1);
    found_tet = *res.begin();
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
    res = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(res.size(), 1);
    found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);

    generate_tri_tet_with_faces(mesh);
//    print_mesh_topology(mesh);
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
    std::set<VertexHandle> expected_tet;
    for (auto expected_tet_vertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expected_tet_vertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[2], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[3]},
            {all_vertices[1], all_vertices[2], all_vertices[3], all_vertices[4]},
    };
    for (auto expected_tet_vertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expected_tet_vertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(res.size(), 5);
    EXPECT_EQ(expected_tets, res);

    /*
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
//    std::set<VertexHandle> expected_tet;
    for (auto expectedTetVertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expectedTetVertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(res.size(), 2);
    EXPECT_EQ(expected_tets, res);
    */

    // Don't detect a tet with some faces missing
    generate_tet_3F(mesh);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_EQ(res.size(), 0);
    res = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(res.size(), 1);
    expected_tets.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[3]}
    };
    expected_tets.clear();
    for (auto expectedTetVertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expectedTetVertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);

    // test some invalid inputs
    mesh.clear();
    res = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(res.empty());
}

// same as test before, but uses find_non_cell_tets_2
TEST_F(TetrahedralMeshBase, find_non_cell_tets_2) {
    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_TRUE(OpenVolumeMesh::find_non_cell_tets_2(mesh, true).empty());
    EXPECT_TRUE(OpenVolumeMesh::find_non_cell_tets_2(mesh, false).empty());

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
    std::set<std::set<VertexHandle>> res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
    EXPECT_TRUE(res.empty());
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(res.size(), 1);
    std::set<VertexHandle> found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);
    std::vector<VertexHandle> face_vertices = {*++mesh.vertices_begin(), *++(++(++mesh.vertices_begin())), *++(++mesh.vertices_begin())};
//    print_mesh_topology(mesh);
    mesh.add_face(face_vertices);
//    print_mesh_topology(mesh);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
    EXPECT_EQ(res.size(), 1);
    found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);

    generate_tet_without_cell(mesh);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
    EXPECT_EQ(res.size(), 1);
    tet_vertices.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        tet_vertices.insert(*v_it);
    }
    found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(res.size(), 1);
    found_tet = *res.begin();
    EXPECT_EQ(tet_vertices, found_tet);

    generate_tri_tet_with_faces(mesh);
//    print_mesh_topology(mesh);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
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
    std::set<VertexHandle> expected_tet;
    for (auto expected_tet_vertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expected_tet_vertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[2], all_vertices[3], all_vertices[4]},
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[3]},
            {all_vertices[1], all_vertices[2], all_vertices[3], all_vertices[4]},
    };
    for (auto expected_tet_vertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expected_tet_vertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(res.size(), 5);
    EXPECT_EQ(expected_tets, res);

    /*
    generate_nested_tets(mesh);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
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
//    std::set<VertexHandle> expected_tet;
    for (auto expectedTetVertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expectedTetVertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(res.size(), 2);
    EXPECT_EQ(expected_tets, res);
     */

    // Don't detect a tet with some faces missing
    generate_tet_3F(mesh);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
    EXPECT_EQ(res.size(), 0);
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(res.size(), 1);
    expected_tets.clear();
    for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        all_vertices.push_back(*v_it);
    }expected_tets_vertices = {
            {all_vertices[0], all_vertices[1], all_vertices[2], all_vertices[3]}
    };
    expected_tets.clear();
    for (auto expectedTetVertices : expected_tets_vertices) {
        expected_tet.clear();
        for (auto expected_vertex : expectedTetVertices) {
            expected_tet.insert(expected_vertex);
        }
        expected_tets.insert(expected_tet);
    }
    EXPECT_EQ(expected_tets, res);

    // test some invalid inputs
    mesh.clear();
    res = OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
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

TEST_F(TetrahedralMeshBase, findNonCellTets_nonFaceTris) {
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
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    // test mixed case
    std::vector<VertexHandle> vertices_cell0(expected_tets.begin()->size());
    std::copy(expected_tets.begin()->begin(), expected_tets.begin()->end(), vertices_cell0.begin());
    mesh.add_cell(vertices_cell0);
    expected_tets.erase(expected_tets.begin());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_EQ(found_tets.size(), 1);
    EXPECT_EQ(found_tets, expected_tets);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

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
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    //check that tet is found even if only_check_faces is true, if the faces exist


    // check that no existing tet is found
    generate_tetrahedral_mesh(mesh);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_TRUE(found_tets.empty());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    // test some invalid inputs
    mesh.clear();
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, false);
    EXPECT_TRUE(found_tets.empty());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());
}


// same as previous test, but uses the second variant of find_non_cell_tets_nonFaceTris
TEST_F(TetrahedralMeshBase, findNonCellTets_2_nonFaceTris) {
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

    auto found_tets = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(found_tets.size(), 2);
    EXPECT_EQ(found_tets, expected_tets);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    // test mixed case
    std::vector<VertexHandle> vertices_cell0(expected_tets.begin()->size());
    std::copy(expected_tets.begin()->begin(), expected_tets.begin()->end(), vertices_cell0.begin());
    mesh.add_cell(vertices_cell0);
    expected_tets.erase(expected_tets.begin());
    found_tets = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(found_tets.size(), 1);
    EXPECT_EQ(found_tets, expected_tets);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

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

    found_tets = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_EQ(found_tets.size(), 3);
    EXPECT_EQ(found_tets, expected_tets);
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    // check that no existing tet is found
    generate_tetrahedral_mesh(mesh);
    found_tets = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_TRUE(found_tets.empty());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());

    // test some invalid inputs
    mesh.clear();
    found_tets = OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    EXPECT_TRUE(found_tets.empty());
    found_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);
    EXPECT_TRUE(found_tets.empty());
}

/*
TEST_F(TetrahedralMeshBase, findNonCellTets_performance) {

    TetrahedralMesh &mesh = this->mesh_;
    generate_tetrahedral_mesh(mesh);

    EXPECT_FALSE(OpenVolumeMesh::contains_void(mesh));

    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.readFile("Cuboid.ovm", mesh);

    auto start = std::clock();
    OpenVolumeMesh::find_non_cell_tets(mesh, true);
    auto end = std::clock();
    auto diff = end - start;
    std::cout << "findNonCellTets time: " << diff << std::endl;
    start = std::clock();
    OpenVolumeMesh::find_non_cell_tets(mesh, false);
    end = std::clock();
    diff = end - start;
    std::cout << "findNonCellTets time: " << diff << std::endl;
    start = std::clock();
    OpenVolumeMesh::find_non_cell_tets_2(mesh, true);
    end = std::clock();
    diff = end - start;
    std::cout << "findNonCellTets_2 time: " << diff << std::endl;
    start = std::clock();
    OpenVolumeMesh::find_non_cell_tets_2(mesh, false);
    end = std::clock();
    diff = end - start;
    std::cout << "findNonCellTets_2 time: " << diff << std::endl;
    EXPECT_TRUE(true);
}
 */

TEST_F(TetrahedralMeshBase, countConnectedComponents) {

    TetrahedralMesh &mesh = this->mesh_;
    mesh.clear();
    EXPECT_EQ(0, OpenVolumeMesh::count_connected_components(mesh));

    generate_tetrahedral_mesh(mesh);
    EXPECT_EQ(1, OpenVolumeMesh::count_connected_components(mesh));


    generate_tets_two_connected_components(mesh);
    EXPECT_EQ(2, OpenVolumeMesh::count_connected_components(mesh));
}

TEST_F(TopologicalLinkBase, faceSetTest_link) {

    TetrahedralMesh tetMesh;
    std::vector<VertexHandle> mesh_vertices;
    generate_triTet(tetMesh, mesh_vertices);

    // manually compute link of one of the vertices incident to all cells
    // first, find one of the two vertices incident to all tets and find the other one, called opposite
    // we can assume that vertices 0 and 4 are incident to all tets, see TopologicalLinkBase::generate_triTet()
    VertexHandle vertex = mesh_vertices[0];
    VertexHandle opposite = mesh_vertices[4];
    // all vertices apart from the original are in the link
    std::set<VertexHandle> vertices = {
            mesh_vertices[1], mesh_vertices[2], mesh_vertices[3], mesh_vertices[4]
    };
    ASSERT_EQ(vertices.size(), 4);
    // the six edges that aren't incident to the original vertex are in the link
    std::set<EdgeHandle> edges;
    for (auto e_it = tetMesh.edges_begin(); e_it != tetMesh.edges_end(); ++e_it) {
        auto heh = tetMesh.halfedge_handle(*e_it, 0);
        if (tetMesh.to_vertex_handle(heh) == vertex ||
            tetMesh.from_vertex_handle(heh) == vertex) {
            continue;
        }
        edges.insert(*e_it);
    }
    ASSERT_EQ(edges.size(), 6);
    // the three faces incident to the opposite vertex are in the link
    std::set<FaceHandle> faces;
    for (auto f_it = tetMesh.faces_begin(); f_it != tetMesh.faces_end(); ++f_it) {
        bool vertex_found = false;
        for (auto fv_it = tetMesh.fv_iter(*f_it); fv_it.valid(); ++fv_it) {
            if ( *fv_it == vertex) {
                vertex_found = true;
                break;
            }
        }
        if (!vertex_found) {
            faces.insert(*f_it);
        }
    }
    ASSERT_EQ(faces.size(), 3);
    std::set<CellHandle> cells;
    // there is no cell in the link
    ASSERT_EQ(cells.size(), 0);

    TopologicalFaceSet expectedLink = {vertices,edges, faces, cells};
    TopologicalFaceSet result = OpenVolumeMesh::link(tetMesh, vertex);
    EXPECT_EQ(expectedLink, result);

    // same for the link of an edge of a tritet, first one of the three edges, which are incident to only one tet
    vertices.clear();
    edges.clear();
    faces.clear();
    cells.clear();
    // choose an edge, which is only incident to one cell, i.e. the valence is 2
    EdgeHandle edge;
    for (auto e_it = tetMesh.edges_begin(); e_it != tetMesh.edges_end(); ++e_it) {
        if (tetMesh.valence(*e_it) == 2) {
            edge = *e_it;
            break;
        }
    }
    auto heh = tetMesh.halfedge_handle(edge, 0);
    // the two vertices not incident to the edge but to the cell are in the link
    auto cell = *tetMesh.ec_iter(edge); // there should only be one incident cell
    for (auto cv_it = tetMesh.cv_iter(cell); cv_it.valid(); ++cv_it) {
        if (*cv_it != tetMesh.to_vertex_handle(heh) &&
            *cv_it != tetMesh.from_vertex_handle(heh)) {
            vertices.insert(*cv_it);
        }
    }
    ASSERT_EQ(vertices.size(), 2);
    // the edge not adjacent to the original edge but incident to the cell is in the link
    cell = *tetMesh.ec_iter(edge); //dito
    for (auto ce_it = tetMesh.ce_iter(cell); ce_it.is_valid(); ++ce_it) {
        auto ce_heh = tetMesh.halfedge_handle(*ce_it, 0);
        if (tetMesh.to_vertex_handle(ce_heh) != tetMesh.to_vertex_handle(heh) &&
            tetMesh.to_vertex_handle(ce_heh)!= tetMesh.from_vertex_handle(heh) &&
            tetMesh.from_vertex_handle(ce_heh) != tetMesh.to_vertex_handle(heh) &&
            tetMesh.from_vertex_handle(ce_heh)!= tetMesh.from_vertex_handle(heh)) {
            edges.insert(*ce_it);
        }
    }
    ASSERT_EQ(edges.size(), 1);
    // There is no face in the link
    ASSERT_EQ(faces.size(), 0);
    // there are no cells in the link
    ASSERT_EQ(cells.size(), 0);

    expectedLink = { vertices, edges, faces, cells};
    result = OpenVolumeMesh::link(tetMesh, edge);
    EXPECT_EQ(expectedLink, result);

    // second, another edge, this time one that is incident to two cells
    vertices.clear();
    edges.clear();
    faces.clear();
    cells.clear();
    // as we know, that the central edge, which is incident to all three cells, is edge(0,4), we can just take edge(0,1)
    bool found = false;
    for (auto out_heh_iter = tetMesh.outgoing_halfedges(mesh_vertices[0]).first; out_heh_iter.is_valid(); ++out_heh_iter) {
        if (tetMesh.to_vertex_handle(*out_heh_iter) == mesh_vertices[1]) {
            edge = tetMesh.edge_handle(*out_heh_iter);
            heh = *out_heh_iter;
            found = true;
            break;
        }
    }
    ASSERT_TRUE(found);
    // all vertices not incident to this edge are in its link
    for (auto v_it = tetMesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        if (*v_it != tetMesh.from_vertex_handle(heh) &&
            *v_it != tetMesh.to_vertex_handle(heh)) {
            vertices.insert(*v_it);
        }
    }
    ASSERT_EQ(vertices.size(), 3);
    // all (both) edges incident only to link vertices and incident to one of the two incident cells are in the link
    // i.e. both edges incident only to link vertices with valence 3
    for (auto e_it = tetMesh.edges_begin(); e_it.is_valid(); ++e_it) {
        auto heh_it = tetMesh.halfedge_handle(*e_it, 0);
        if (((vertices.find(tetMesh.to_vertex_handle(heh_it))) != vertices.end()) &&
            ((vertices.find(tetMesh.from_vertex_handle(heh_it))) != vertices.end()) &&
            (tetMesh.valence(*e_it) == 3)) {
            edges.insert(*e_it);
        }
    }
    ASSERT_EQ(edges.size(), 2);
    // There are no faces in the link
    ASSERT_EQ(faces.size(), 0);
    // There are no cells in the link
    ASSERT_EQ(cells.size(), 0);
    expectedLink = {vertices, edges, faces, cells};
    result = link(tetMesh, edge);
    EXPECT_EQ(expectedLink, result);

    // third, even another edge, this time the one with three incident cells
    vertices.clear();
    edges.clear();
    faces.clear();
    cells.clear();
    // again, we know that the edge(0,4) is the one we are looking for
    found = false;
    for (auto out_heh_iter = tetMesh.outgoing_halfedges(mesh_vertices[0]).first; out_heh_iter.is_valid(); ++out_heh_iter) {
        if (tetMesh.to_vertex_handle(*out_heh_iter) == mesh_vertices[4]) {
            edge = tetMesh.edge_handle(*out_heh_iter);
            heh = *out_heh_iter;
            found = true;
            break;
        }
    }
    ASSERT_TRUE(found);
    // all vertices not incident to the edge are in the link
    for (auto v_it = tetMesh.vertices_begin(); v_it.is_valid(); ++v_it) {
        if (*v_it != tetMesh.from_vertex_handle(heh) &&
            *v_it != tetMesh.to_vertex_handle(heh)) {
            vertices.insert(*v_it);
        }
    }
    ASSERT_EQ(vertices.size(), 3);
    // all edges not adjacent to the edge are in the link
    for (auto e_it = tetMesh.edges_begin(); e_it.is_valid(); ++e_it) {
        auto heh_it = tetMesh.halfedge_handle(*e_it, 0);
        if (tetMesh.from_vertex_handle(heh) != tetMesh.from_vertex_handle(heh_it) &&
            tetMesh.from_vertex_handle(heh) != tetMesh.to_vertex_handle(heh_it) &&
            tetMesh.to_vertex_handle(heh) != tetMesh.from_vertex_handle(heh_it) &&
            tetMesh.to_vertex_handle(heh) != tetMesh.to_vertex_handle(heh_it)) {
            edges.insert(*e_it);
        }
    }
    ASSERT_EQ(edges.size(), 3);
    // there are no faces in the link
    ASSERT_EQ(faces.size(), 0);
    // there are no cells in the link
    ASSERT_EQ(cells.size(), 0);
    expectedLink = { vertices, edges, faces, cells};
    result = link(tetMesh, edge);
    EXPECT_EQ(expectedLink, result);
}