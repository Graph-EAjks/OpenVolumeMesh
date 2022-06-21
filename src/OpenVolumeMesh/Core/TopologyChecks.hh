#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Core/TopologicalLink.hh>

namespace OpenVolumeMesh{

    CellHandle cell_exists(const TetrahedralMeshTopologyKernel& mesh,
                                  const std::vector<VertexHandle>& cell_vertices);

    bool face_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const FaceHandle& face);

    bool cell_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const CellHandle& cell);

    // TODO: This method does not work in all cases. It only finds non-cell tets which are completely surrounded by faces.
    std::set<std::set<VertexHandle>> findNonCellTets(TetrahedralMeshTopologyKernel& mesh);

    bool link_condition(const TetrahedralMeshTopologyKernel& mesh,
                        const EdgeHandle& edge);

    bool singleConnectedComponent(TetrahedralMeshTopologyKernel&  mesh);

    bool containsVoid(TetrahedralMeshTopologyKernel&  mesh);

    std::vector<VertexHandle> nonManifoldBoundaryVertices(TetrahedralMeshTopologyKernel& mesh);

    bool manifoldVertex(TetrahedralMeshTopologyKernel& mesh,
                        const VertexHandle& vertex);

    bool noDoubleEdges(TetrahedralMeshTopologyKernel& mesh);

    void printMeshTopology(TetrahedralMeshTopologyKernel& mesh);
}