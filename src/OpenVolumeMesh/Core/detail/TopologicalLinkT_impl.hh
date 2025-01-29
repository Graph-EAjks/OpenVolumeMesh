#include "TopologicalLink.hh"


namespace OpenVolumeMesh {

    template<class MeshT>
    TopologicalFaceSet link(const MeshT &mesh,
                            const VertexHandle vertex) {

        VertexSet vertices;
        EdgeSet edges;
        FaceSet faces;

        // all neighbouring vertices are in the link
        for (auto v = mesh.vv_iter(vertex); v.valid(); v++) {
            vertices.insert(*v);
        }

        // all edges, which are opposite of the vertex in a face are in the link and all faces which are opposite of
        // the vertex in a tet including all their edges are in the link.
        for (auto c = mesh.vc_iter(vertex); c.valid(); c++) {

            for (auto hf = mesh.chf_iter(*c); hf.valid(); hf++) {
                bool found(false);
                for (auto v = mesh.hfv_iter(*hf); v.valid(); v++) {
                    if (*v == vertex) {
                        found = true;
                        break;
                    }
                }

                // the face does not contain the original vertex, so it is on the opposite side of the tet
                if (!found) {
                    faces.insert(mesh.face_handle(*hf));
                    for (auto e = mesh.hfe_iter(*hf); e.valid(); e++) {
                        edges.insert(*e); // an edge might be inserted more than once, but that's ok as edges is a set
                    }
                }
            }
        }

        // also include edges which are not part of an incident cell of the vertex, but only of an incident face
        for (auto f = mesh.vf_iter(vertex); f.valid(); f++) {
            for (auto e = mesh.fe_iter(*f); e.valid(); e++) {
                bool found(false);
                if (mesh.edge(*e).to_vertex() == vertex ||
                    mesh.edge(*e).from_vertex() == vertex) {
                    found = true;
                }
                if (!found) {
                    edges.insert(*e);
                }
            }
        }

        return {vertices, edges, faces, {}};
    }


    template<class MeshT>
    TopologicalFaceSet link(const MeshT &mesh,
                            const EdgeHandle edge) {

        VertexSet vertices;
        EdgeSet edges;

        auto from_vertex = mesh.edge(edge).from_vertex();
        auto to_vertex = mesh.edge(edge).to_vertex();

        // all edges opposite of a tet including their end vertices are in the link
        for (auto c = mesh.ec_iter(edge); c.valid(); c++) {

            for (auto e = mesh.ce_iter(*c); e.valid(); e++) {

                auto other_from_vertex = mesh.edge(*e).from_vertex();
                auto other_to_vertex = mesh.edge(*e).to_vertex();

                // check if found edge is adjacent to the original edge
                if (from_vertex != other_from_vertex &&
                    from_vertex != other_to_vertex &&
                    to_vertex != other_from_vertex &&
                    to_vertex != other_to_vertex) {

                    edges.insert(*e);
                    vertices.insert(other_from_vertex);
                    vertices.insert(other_to_vertex);
                }
            }
        }

        // also check for neighbouring vertices not incident to the same cell
        for (auto f = mesh.ef_iter(edge); f.valid(); f++) {
            for (auto v = mesh.fv_iter(*f); v++;) {
                if (*v != mesh.edge(edge).to_vertex() &&
                    *v != mesh.edge(edge).from_vertex()) {
                    vertices.insert(*v);
                }
            }
        }

        return {vertices, edges, {}, {}};

    }


    template<class MeshT>
    TopologicalFaceSet link_outsiders(const MeshT &mesh,
                                      const OpenVolumeMesh::EdgeHandle edge) {

        auto edge_link = link(mesh, edge);
        auto from_vertex_link = link(mesh, mesh.edge(edge).from_vertex());
        auto to_vertex_link = link(mesh, mesh.edge(edge).to_vertex());

        return from_vertex_link.intersection(to_vertex_link).subtract(edge_link);
    }


    template<class MeshT>
    bool link_condition(const MeshT &mesh,
                        const EdgeHandle edge) {

        auto from_vertex = mesh.edge(edge).from_vertex();
        auto to_vertex = mesh.edge(edge).to_vertex();

        auto vertex_link_intersection = link(mesh, from_vertex).intersection(link(mesh, to_vertex));

        //the link of an edge cannot contain any face.
        // Therefore if the intersection's link contains a face, it cannot be equal to the edge's link
        if (vertex_link_intersection.faces().size()) {
            return false;
        }

        return vertex_link_intersection == link(mesh, edge);
    }


    template<class MeshT>
    bool link_condition(const MeshT &mesh,
                        const OpenVolumeMesh::HalfEdgeHandle heh) {

        return link_condition(mesh, mesh.edge_handle(heh));
    }

}