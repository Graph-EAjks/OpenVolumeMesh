#include "Topo_checks_examples.hh"
#include "OpenVolumeMesh/FileManager/FileManager.hh"

int main(int _argc, char** _argv) {

    bool success = true;

    success &= test_cell_exists();
    success &= test_face_contains_vertex();
    success &= test_cell_contains_vertex();
    success &= test_find_non_cell_tets();
    success &= test_single_connected_component();
    success &= test_count_connected_components();
    success &= test_contains_void();

    if (success) {
        std::cout << "success!" << std::endl;
    } else {
        std::cout << "something went wrong" << std::endl;
    }

    return 0;
}

bool test_cell_exists() {
    bool success = true;
    // We generate a simple tet mesh consisting of one cell
    generate_simple_tetmesh(mesh);

    // All vertices should define exactly one tet, which should be a cell
    auto cell_vertices = mesh.get_cell_vertices(*mesh.cells_begin());

    // There is only one cell in this mesh which is the cell that we try to find
    auto cell_expected = *begin(mesh.cells());

    // returns a valid cell, if cell_vertices define a tet which wis a cell
    auto cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    success &= cell.has_value();
    success &= cell == cell_expected;

    // Now, we generate another mesh consisting of one tet, which is not a cell, so no cell should be found
    generate_simple_tetmesh_without_cell(mesh);

    // There are four vertices in this mesh which form a tet. So, one could expect that they form a cell.
    cell_vertices.clear();
    for (auto vertex : mesh.vertices()) {
        cell_vertices.push_back(vertex);
    }

    // However, the cell does not exist.
    cell = OpenVolumeMesh::cell_exists(mesh, cell_vertices);
    success &= !cell.has_value();

    if (!success) {
        std::cout << "test_cell_exists failed" << std::endl;
    }
    return success;
}

bool test_face_contains_vertex() {
    bool success = true;

    // We generate a tetmesh containing one face
    generate_triangle(mesh);

    // There is exactly one face in this mesh
    auto face = *mesh.faces_begin();

    // And this face contains all vertices in the mesh
    auto vertex = *mesh.vertices_begin();
    success &= OpenVolumeMesh::face_contains_vertex(mesh, vertex, face);

    // Now, let's add another vertex which is not in this face
    VertexHandle other_vertex = mesh.add_vertex();

    // The new vertex is not in the face
    success &= !OpenVolumeMesh::face_contains_vertex(mesh, other_vertex, face);

    if (!success) {
        std::cout << "test_face_contains_vertex failed" << std::endl;
    }
    return success;
}

bool test_cell_contains_vertex() {
    bool success = true;

    // A simple mesh with one tet
    generate_simple_tetmesh(mesh);

    // There is only one cell in this mesh
    auto cell = *mesh.cells_begin();

    // Every vertex is in the only cell of the mesh
    auto vertex = *mesh.vertices_begin();
    success &= OpenVolumeMesh::cell_contains_vertex(mesh, vertex, cell);

    // Now, we add another vertex which is not part of the cell
    auto other_vertex = mesh.add_vertex();
    success &= !OpenVolumeMesh::cell_contains_vertex(mesh, other_vertex, cell);

    if (!success) {
        std::cout << "test_cell_contains_vertex failed" << std::endl;
    }
    return success;
}

bool test_find_non_cell_tets() {
    bool success = true;

    // A simple mesh with one tet, which is a cell
    generate_simple_tetmesh(mesh);

    auto non_cell_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);

    //There is only one tet, which is a cell, so there should be no non-cell-tets
    success &= non_cell_tets.empty();

    // Another simple mesh with one tet, which is not a cell
    generate_simple_tetmesh_without_cell(mesh);

    non_cell_tets = OpenVolumeMesh::find_non_cell_tets(mesh, true);

    // Now, the tet, which is not a cell should be found
    success &= !non_cell_tets.empty();

    // Check, that it found exactly one tet with four vertices
    success &= non_cell_tets.size() == 1;
    auto tet = *non_cell_tets.begin();
    success &= tet.size() == 4;

    if (!success) {
        std::cout << "test_find_non_cell_tets failed" << std::endl;
    }
    return success;
}

bool test_single_connected_component() {
    bool success = true;

    // An empty mesh has not exactly one connected component, so it should return false
    generate_empty_mesh(mesh);
    success &= !OpenVolumeMesh::single_connected_component(mesh);

    // A mesh containing exactly one tet is connected
    generate_simple_tetmesh(mesh);
    success &= OpenVolumeMesh::single_connected_component(mesh);

    // A mesh containing two unconnected tets is not connected
    generate_two_unconnected_tets(mesh);
    success &= !OpenVolumeMesh::single_connected_component(mesh);

    if (!success) {
        std::cout << "test_single_connected_component failed" << std::endl;
    }

    return success;
}

bool test_count_connected_components() {
    bool success = true;

    // An empty mesh has 0 connected components
    generate_empty_mesh(mesh);
    success &= OpenVolumeMesh::count_connected_components(mesh) == 0;

    // A mesh containing exactly 1 tet has 1 connected components
    generate_simple_tetmesh(mesh);
    success &= OpenVolumeMesh::count_connected_components(mesh) == 1;

    // A mesh containing 2 tets has 2 conenctec components
    generate_two_unconnected_tets(mesh);
    success &= OpenVolumeMesh::count_connected_components(mesh) == 2;

    if (!success) {
        std::cout << "test_count_connected_components failed" << std::endl;
    }

    return success;
}

bool test_contains_void() {
    bool success = true;

    // A simple tet should not contain any void
    generate_simple_tetmesh(mesh);
    success &= !OpenVolumeMesh::contains_void(mesh);

    generate_mesh_with_void(mesh);
    success &= OpenVolumeMesh::contains_void(mesh);

    if (!success) {
        std::cout << "test_contains_void failed" << std::endl;
    }
    return success;
}

/**
 * We now define some mesh generators. As we are only interested in the topology and not in the
 * geometry, we don't need to specify coordinates for the vertices.
 */

void generate_empty_mesh(TetrahedralMesh& _mesh) {
    _mesh.clear();
}

void generate_triangle(TetrahedralMesh& _mesh) {
    _mesh.clear();

    // Create vertices
    VertexHandle v1 = _mesh.add_vertex();
    VertexHandle v2 = _mesh.add_vertex();
    VertexHandle v3 = _mesh.add_vertex();

    // Create face
    std::vector<VertexHandle> face_vertices = {v1, v2, v3};
    _mesh.add_face(face_vertices);
}

/**
 * Generates a simple tetrahedral mesh consisting of one tet.
 */
void generate_simple_tetmesh(TetrahedralMesh& _mesh) {

    _mesh.clear();

    // Create vertices
    VertexHandle v1 = _mesh.add_vertex();
    VertexHandle v2 = _mesh.add_vertex();
    VertexHandle v3 = _mesh.add_vertex();
    VertexHandle v4 = _mesh.add_vertex();

    // Add  cell
    _mesh.add_cell(v1, v2, v3, v4);
}

/**
 * Generates a simple tetrahedral mesh consisting of one tet, which is not a cell
 */
void generate_simple_tetmesh_without_cell(TetrahedralMesh& _mesh){

    _mesh.clear();

    // Create vertices
    VertexHandle v1 = _mesh.add_vertex();
    VertexHandle v2 = _mesh.add_vertex();
    VertexHandle v3 = _mesh.add_vertex();
    VertexHandle v4 = _mesh.add_vertex();

    // Add faces. This also creates edges, halfedges and halffaces, but no cell
    // The order of vertices in the faces is chosen like this to ensure a correct orientation.
    _mesh.add_face({v1, v2, v3});
    _mesh.add_face({v1, v4, v2});
    _mesh.add_face({v1, v3, v4});
    _mesh.add_face({v2, v4, v3});
}

void generate_two_unconnected_tets(TetrahedralMesh& _mesh) {

    _mesh.clear();

    // Create vertices
    VertexHandle vertices[8];
    for (int i = 0; i < 8; ++i) {
        vertices[i] = _mesh.add_vertex();
    }

    // Create cells
    _mesh.add_cell(vertices[0], vertices[1], vertices[2], vertices[3]);
    _mesh.add_cell(vertices[4], vertices[5], vertices[6], vertices[7]);
}

void generate_mesh_with_void(TetrahedralMesh& _mesh) {

    _mesh.clear();

    // First, we import a large mesh which contains several non-boundary cells
    OpenVolumeMesh::IO::FileManager fileManager;
    //TODO: "End of file reached while searching for input!" is printed, but it seems it works correctly
    fileManager.readFile("Files/Cuboid.ovm", _mesh);

    // Then, we look for a non-boundary cell and remove it, so that a void is created
    CellHandle cell_to_remove(-1);
    for (auto cell : _mesh.cells()){
        if (_mesh.is_boundary(cell)) continue;
        bool has_boundary_vertex = false;
        for (auto cv_it = _mesh.cv_iter(cell); cv_it.is_valid(); ++cv_it) {
            if (_mesh.is_boundary(*cv_it)) {
                has_boundary_vertex = true;
                break;
            }
        }
        if (!has_boundary_vertex) {
            cell_to_remove = cell;
            break;
        }
    }

    if(!cell_to_remove.is_valid()) {
        _mesh.clear();
        return;
    }

    _mesh.delete_cell(cell_to_remove);
}