#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Core/TopologicalLink.hh>
#include <OpenVolumeMesh/Core/detail/exceptions.hh>

namespace OpenVolumeMesh{

    //TODO: check, which functions rely on tet structure and which can be used for any mesh

    CellHandle cell_exists(const TetrahedralMeshTopologyKernel& mesh,
                                  const std::vector<VertexHandle>& cell_vertices);

    bool face_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const FaceHandle& face);

    bool cell_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const CellHandle& cell);

    /**
     *
     * @param mesh
     * @param only_check_faces only returns non_cell_tets, where the four faces exist, otherwise also return tets with
     * missing faces
     * @return
     */
    std::set<std::set<VertexHandle>> find_non_cell_tets(const TetrahedralMeshTopologyKernel& mesh, bool only_check_faces);

    std::set<std::set<VertexHandle>> find_non_cell_tets_2(const TetrahedralMeshTopologyKernel& mesh, bool only_check_faces);

    bool link_condition(const TetrahedralMeshTopologyKernel& mesh,
                        const EdgeHandle& edge);

    bool single_connected_component(const TetrahedralMeshTopologyKernel&  mesh);

    size_t count_connected_components(const TetrahedralMeshTopologyKernel& mesh);

    bool contains_void(const TetrahedralMeshTopologyKernel&  mesh);

    bool manifold_vertex(const TetrahedralMeshTopologyKernel& mesh,
                         const VertexHandle& vertex);

    bool no_double_edges(const TetrahedralMeshTopologyKernel& mesh);

    std::set<std::set<HEH>> find_multi_edges(const TetrahedralMeshTopologyKernel& mesh);

    void print_mesh_topology(const TetrahedralMeshTopologyKernel& mesh);

    std::set<std::set<VertexHandle>> find_non_face_triangles(const TetrahedralMeshTopologyKernel& mesh);

    std::set<std::set<VertexHandle>> find_non_face_triangles_around_vertex(const TetrahedralMeshTopologyKernel& mesh, const VertexHandle& vertex);
}