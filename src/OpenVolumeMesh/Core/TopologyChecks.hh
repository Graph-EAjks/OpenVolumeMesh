#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>

namespace OpenVolumeMesh{

    OVM_EXPORT
    CellHandle cell_exists(const TetrahedralMeshTopologyKernel& mesh,
                                  const std::vector<VertexHandle>& cell_vertices);


    OVM_EXPORT
    bool face_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const FaceHandle& face);

    OVM_EXPORT
    bool cell_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                     const VertexHandle& vertex,
                                     const CellHandle& cell);

    /*
    OVM_EXPORT
    std::set<std::pair<std::set<VertexHandle>, bool>> findNonCellTets(const TetrahedralMeshTopologyKernel&);

    OVM_EXPORT
    bool link_condition(const TetrahedralMeshTopologyKernel& mesh,
                        const EdgeHandle& edge);
    */

    bool singleConnectedComponent(TetrahedralMeshTopologyKernel&  mesh);

    bool containsVoid(TetrahedralMeshTopologyKernel&  mesh);

    std::vector<VertexHandle> nonManifoldBoundaryVertices(TetrahedralMeshTopologyKernel& mesh);

    bool manifoldVertex(TetrahedralMeshTopologyKernel& mesh,
                        const VertexHandle& vertex);

    bool noDoubleEdges(TetrahedralMeshTopologyKernel& mesh);
}