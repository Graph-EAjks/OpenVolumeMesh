#ifndef OPENVOLUMEMESH_TOPO_CHECKS_EXAMPLES_HH
#define OPENVOLUMEMESH_TOPO_CHECKS_EXAMPLES_HH

// C++ includes
#include <iostream>
#include <vector>

#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/TopologyChecks.hh>

typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d TetrahedralMesh;
typedef OpenVolumeMesh::VertexHandle VertexHandle;

using namespace OpenVolumeMesh;

bool test_cell_exists();
bool test_face_contains_vertex();
bool test_cell_contains_vertex();
bool test_find_non_cell_tets();
bool test_single_connected_component();
bool test_count_connected_components();
bool test_contains_void();

void generate_empty_mesh(TetrahedralMesh& _mesh);
void generate_triangle(TetrahedralMesh& _mesh);
void generate_simple_tetmesh(TetrahedralMesh& _mesh);
void generate_simple_tetmesh_without_cell(TetrahedralMesh& _mesh);
void generate_two_unconnected_tets(TetrahedralMesh& _mesh);
void generate_mesh_with_void(TetrahedralMesh& _mesh);

TetrahedralMesh mesh;

#endif //OPENVOLUMEMESH_TOPO_CHECKS_EXAMPLES_HH
