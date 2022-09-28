#pragma once

#include "TopologicalFaceSet.hh"

/**
 * This class offers functions in connection with the topological link condition. This concept is described in detail in
 * https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.40.7253
 *
 * It is basically used to determine, if an edge in a tetrahedral mesh can be collapsed without making the mesh
 * topologically degenerate. The link condition for an edge is a necessary and sufficient condition to be able to
 * collapse it without problems.
 *
 * The class TopologicalFaceSet is only used by this class in order to facilitate handling links (a set of vertices,
 * edges, faces and cells).
 */

namespace OpenVolumeMesh {
/* Computes the 'Link' of a vertex in the topological sense */
    template<class MeshT>
    TopologicalFaceSet link(const MeshT &mesh,
                            const OpenVolumeMesh::VertexHandle vertex);


/* Computes the 'Link' of an edge in the topological sense */
    template<class MeshT>
    TopologicalFaceSet link(const MeshT &mesh,
                            const OpenVolumeMesh::EdgeHandle edge);


/* Returns the set of topological faces that are in the intersection of
the edge's vertices' links but not in the edge's link */
    template<class MeshT>
    TopologicalFaceSet link_outsiders(const MeshT &mesh,
                                      const OpenVolumeMesh::EdgeHandle edge);


/* Checks that the 'link condition' for a particular edge holds.
i.e. If edge goes from vertex a to vertex b, it checks whether

Lk(a) [intersection] Lk(b) = Lk(ab)

is true or not.

This can typically be used to check if this edge can be safely collapsed
without hurting the mesh's topology. */
    template<class MeshT>
    bool link_condition(const MeshT &mesh,
                        const OpenVolumeMesh::EdgeHandle eh);


    template<class MeshT>
    bool link_condition(const MeshT &mesh,
                        const OpenVolumeMesh::HalfEdgeHandle heh);
}