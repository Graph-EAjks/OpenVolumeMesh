#include <OpenVolumeMesh/Core/TopologyChecks.hh>
#include <OpenVolumeMesh/Core/TopologicalLinkT_impl.hh>

namespace OpenVolumeMesh {

    CellHandle cell_exists(const TetrahedralMeshTopologyKernel& mesh,
                           const std::vector<VertexHandle>& cell_vertices) {
        if(cell_vertices.size() == 4){
            for(auto vc_it = mesh.vc_iter(cell_vertices[0]); vc_it.valid(); vc_it++){
                int found_vertices_count(0);
                for(auto cv_it = mesh.cv_iter(*vc_it); cv_it.valid(); cv_it++){
                    if(std::find(cell_vertices.begin(), cell_vertices.end(), *cv_it) != cell_vertices.end()){
                        found_vertices_count++;
                    }
                }

                if(found_vertices_count == 4){
                    return *vc_it;
                }
            }
        }
        return CellHandle(-1);
    }

    bool face_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                          const VertexHandle& vertex,
                                          const FaceHandle& face){
        if (!mesh.is_valid(vertex) || !mesh.is_valid(face)) {
            return false;
        }
        auto hf_vertices = mesh.get_halfface_vertices(mesh.halfface_handle(face, 0));

        return  hf_vertices[0] == vertex ||
                hf_vertices[1] == vertex ||
                hf_vertices[2] == vertex;
    }

    bool cell_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                          const VertexHandle& vertex,
                                          const CellHandle& cell){
        if (!mesh.is_valid(vertex) || !mesh.is_valid(cell)) {
            return false;
        }
        auto c_vertices = mesh.get_cell_vertices(cell);

        return  c_vertices[0] == vertex ||
                c_vertices[1] == vertex ||
                c_vertices[2] == vertex ||
                c_vertices[3] == vertex;
    }

    std::set<std::pair<std::set<VertexHandle>, bool>> findNonCellTets(TetrahedralMeshTopologyKernel& mesh_){

        std::set<std::pair<std::set<VertexHandle>, bool>> non_cell_tets;

        for(auto v: mesh_.vertices()){

            auto neighbor = mesh_.request_vertex_property<bool>("neighbor");

            for(auto vv_it = mesh_.vv_iter(v); vv_it.valid(); vv_it++){
                neighbor[*vv_it] = true;
            }

            for(auto vf_it = mesh_.vf_iter(v); vf_it.valid(); vf_it++){

                std::set<VertexHandle> tet_vertices;
                tet_vertices.insert(v);

                //find opposite edge to v
                EdgeHandle opposite_edge(-1);

                for(auto fe_it = mesh_.fe_iter(*vf_it); fe_it.valid(); fe_it++){
                    if(mesh_.edge(*fe_it).from_vertex() != v &&
                       mesh_.edge(*fe_it).to_vertex() != v){
                        tet_vertices.insert(mesh_.edge(*fe_it).from_vertex());
                        tet_vertices.insert(mesh_.edge(*fe_it).to_vertex());

                        opposite_edge = *fe_it;
                    }
                }

                //iterate through the opposite edge's faces to try to find another neighbor vertex
                for(auto ef_it = mesh_.ef_iter(opposite_edge); ef_it.valid(); ef_it++){
                    std::set<VertexHandle> partial_cell = tet_vertices;

                    //find the opposite vertex
                    for(auto fv_it = mesh_.fv_iter(*ef_it); fv_it.valid(); fv_it++){
                        if(mesh_.edge(opposite_edge).from_vertex() != *fv_it &&
                           mesh_.edge(opposite_edge).to_vertex() != *fv_it &&
                           neighbor[*fv_it]){

                            partial_cell.insert(*fv_it);

                            bool found_cell(false);

                            //check if it's a cell or not
                            for(auto vc_it = mesh_.vc_iter(v); vc_it.valid(); vc_it++){
                                int found_vertices(0);
                                for(auto cv_it = mesh_.cv_iter(*vc_it); cv_it.valid(); cv_it++){
                                    if(partial_cell.find(*cv_it) != partial_cell.end()){
                                        found_vertices++;
                                    }
                                }
                                if(found_vertices == 4){
                                    found_cell = true;
                                }
                            }

                            if(!found_cell){

                                if(partial_cell.size() != 4){
                                    std::cout<<" ERROR - tet does not contain 4 vertices"<<std::endl;
                                    return {};
                                }
                                //check that at least two vertices are non-boundary
                                int boundary_count(0);
                                for(auto v: partial_cell){
                                    boundary_count += mesh_.is_boundary(v);
                                }
                                if(boundary_count <=2){

                                    bool all_edge_non_collapsible(true);
                                    std::set<EdgeHandle> tet_edges;
                                    for(auto v: partial_cell){
                                        for(auto out_he: mesh_.outgoing_halfedges(v)){
                                            auto edge = mesh_.edge_handle(out_he);

                                            if(partial_cell.find(mesh_.to_vertex_handle(out_he)) != partial_cell.end()){
                                                tet_edges.insert(edge);
                                                all_edge_non_collapsible &= !link_condition(mesh_, edge);
                                            }
                                        }
                                    }

                                    if(tet_edges.size() != 6){
                                        std::cout<<" ERROR - tet does not contain 6 edges"<<std::endl;
                                        return {};
                                    }

                                    std::set<FaceHandle> tet_faces;
                                    for(auto e: tet_edges){
                                        for(auto ef_it = mesh_.ef_iter(e); ef_it.valid(); ef_it++){
                                            for(auto fv_it = mesh_.fv_iter(*ef_it); fv_it.valid(); fv_it++){
                                                if(*fv_it != mesh_.edge(e).from_vertex() &&
                                                   *fv_it != mesh_.edge(e).to_vertex() &&
                                                   partial_cell.find(*fv_it) != partial_cell.end()){
                                                    tet_faces.insert(*ef_it);
                                                }
                                            }
                                        }
                                    }

                                    if(tet_faces.size() == 4){
                                        non_cell_tets.insert({partial_cell, all_edge_non_collapsible});
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return non_cell_tets;
    }


    bool link_condition(const TetrahedralMeshTopologyKernel& mesh,
                        const EdgeHandle& edge){

        auto from_vertex = mesh.edge(edge).from_vertex();
        auto to_vertex = mesh.edge(edge).to_vertex();

        auto vertex_link_intersection = link(mesh, from_vertex).intersection(link(mesh, to_vertex));

        if(vertex_link_intersection.faces().size()){
            return false;
        }

        return vertex_link_intersection == link(mesh, edge);
    }

    bool singleConnectedComponent(TetrahedralMeshTopologyKernel&  mesh){

        if (!mesh.vertices_begin().is_valid()) { // mesh contains no vertices, so it is empty
            return false; //TODO: it might be preferable to return true, as an empty mesh is connected (maybe not?)
        }

        auto visited_prop = mesh.request_vertex_property<bool>("visited");

        std::set<VertexHandle> to_visit;
        to_visit.insert(*mesh.vertices().first);

        while(!to_visit.empty()){
            auto next = *to_visit.begin();
            to_visit.erase(to_visit.begin());

            visited_prop[next] = true;

            for(auto out_he_it = mesh.outgoing_halfedges(next).first; out_he_it.valid(); out_he_it++){
                auto to_vertex = mesh.to_vertex_handle(*out_he_it);
                if(!visited_prop[to_vertex]){
                    to_visit.insert(to_vertex);
                }
            }
        }

        bool all_visited(true);

        for(auto v: mesh.vertices()){
            if(!visited_prop[v]){
                all_visited = false;
            }
        }


        return all_visited;
    }

    bool containsVoid(TetrahedralMeshTopologyKernel&  mesh){

        std::vector<FaceHandle> to_visit;

        for(auto f: mesh.faces()){
            if(mesh.is_boundary(f)){
                to_visit.push_back(f);
                break;
            }
        }

        if(!to_visit.size()){
            std::cerr<<" couldn't find a boundary face, which is very weird"<<std::endl;
            return false;
        }

        auto visited = mesh.request_face_property<bool>("visited");

        int iteration_count(0);

        while(to_visit.size()){
            //std::cout<<" stack size at iteration "<<iteration_count<<": "<<to_visit.size()<<std::endl;

            auto current_face = to_visit.back();
            to_visit.pop_back();

            if(!visited[current_face]){

                visited[current_face] = true;

                for(auto fe_it = mesh.fe_iter(current_face); fe_it.valid(); fe_it++){

                    for(auto ef_it = mesh.ef_iter(*fe_it); ef_it.valid(); ef_it++){
                        if(mesh.is_boundary(*ef_it) && !visited[*ef_it]){
                            to_visit.push_back(*ef_it);
                        }
                    }
                }
            }

            iteration_count++;
        }


        for(auto f: mesh.faces()){
            if(mesh.is_boundary(f) && !visited[f]){
                std::cout<<" found a face not visited"<<std::endl;
                return true;
            }
        }

        return false;
    }

    std::vector<VertexHandle> nonManifoldBoundaryVertices(TetrahedralMeshTopologyKernel& mesh){

        std::vector<VertexHandle> non_manifold_boundary_vertices;

        for(auto v: mesh.vertices()){
            if(mesh.is_boundary(v)){

                if(!manifoldVertex(mesh, v)){
                    non_manifold_boundary_vertices.push_back(v);
                }

            }
        }

        return non_manifold_boundary_vertices;
    }

    bool manifoldVertex(TetrahedralMeshTopologyKernel& mesh,
                                    const VertexHandle& vertex){


        if (!mesh.is_valid(vertex)) {
            return false;
        }

        // if the vertex is incident to less than three edges, it cannot be part of a tet
        if (mesh.valence(vertex) < 3) {
            return false;
        }
        // if an incident edge is incident to less than two faces, it cannot be part of a tet
        for(auto out_he: mesh.outgoing_halfedges(vertex)){
            auto eh = mesh.edge_handle(out_he);
            if (mesh.valence(eh) < 2) {
                return false;
            }
        }
        // if an incident face is incident to no cell, then it is not part of a tet
        for (auto vc_it = mesh.vc_iter(vertex); vc_it.is_valid(); ++vc_it) {
            if (mesh.valence(*vc_it) < 1) {
                return false;
            }
        }

        if(!mesh.is_boundary(vertex)){
            return true;
        }

        //find an initial boundary face
        FaceHandle initial_boundary_face(-1);
        for(auto vf_it = mesh.vf_iter(vertex); vf_it.valid(); vf_it++){
            if(mesh.is_boundary(*vf_it)){
                initial_boundary_face = *vf_it;
                break;
            }
        }

        if(initial_boundary_face.idx() == -1){
            std::cerr<<" couldn't find boundary face adjacent to boundary vertex. returning false"<<std::endl;
            return false;
        }

        //std::cout<<" ----------------------------"<<std::endl;
        //std::cout<<" checking vertex "<<vertex<<std::endl;

        //std::cout<<" -- initial boundary face : "<<initial_boundary_face<<": "<<mesh.face(initial_boundary_face)<<std::endl;


        auto visited_edge = mesh.request_edge_property<bool>("visited");
        auto visited_face = mesh.request_face_property<bool>("visited");

        EdgeHandle next_edge(-1);

        //visit all vertices of the initial face
        for(auto fe_it = mesh.fe_iter(initial_boundary_face); fe_it.valid(); fe_it++){
            if(mesh.edge(*fe_it).to_vertex() == vertex ||
               mesh.edge(*fe_it).from_vertex() == vertex){
                visited_edge[*fe_it] = true;
                //std::cout<<" ---- visited edge "<<*fe_it<<" and marked it as current edge"<<std::endl;
                next_edge = *fe_it;
            }
        }
        visited_face[initial_boundary_face] = true;

        while(next_edge.idx() != -1){
            //std::cout<<" --------- "<<std::endl;

            auto current_edge = next_edge;
            next_edge = EdgeHandle(-1);

            //circulate through boundary faces around edge to find the next one
            for(auto ef_it = mesh.ef_iter(current_edge); ef_it.valid(); ef_it++){
                //std::cout<<" -- checking face "<<*ef_it<<": "<<mesh.face(*ef_it)<<std::endl;
                if(mesh.is_boundary(*ef_it) && ! visited_face[*ef_it]){
                    //found a neighboring face to both current_edge and vertex
                    //so visit the next edge and use it as next current_edge

                    //std::cout<<" ---- neighbor face: "<<*ef_it<<": "<<mesh.face(*ef_it)<<std::endl;

                    for(auto fe_it = mesh.fe_iter(*ef_it); fe_it.valid(); fe_it++){
                        if(mesh.is_boundary(*fe_it)&&
                           (mesh.edge(*fe_it).to_vertex() == vertex ||
                            mesh.edge(*fe_it).from_vertex() == vertex) &&
                           !visited_edge[*fe_it]){
                            visited_edge[*fe_it] = true;
                            next_edge = *fe_it;
                            //std::cout<<" ---- next neighbor edge: "<<*fe_it<<": "<<mesh.edge(*fe_it)<<std::endl;
                        }
                    }
                }
                visited_face[*ef_it] = true;
            }
        }

        //and finally we check that all boundary faces were visited
        bool all_visited(true);
        for(auto vf_it = mesh.vf_iter(vertex); vf_it.valid(); vf_it++){
            if(mesh.is_boundary(*vf_it)){
                all_visited &= visited_face[*vf_it];
            }
        }
        if(!all_visited){
            std::cout<<" vertex "<<vertex<<" is non-manifold"<<std::endl;
        }

        return all_visited;
    }

    bool noDoubleEdges(TetrahedralMeshTopologyKernel& mesh){

        for(auto v: mesh.vertices()){
            auto visited = mesh.request_vertex_property<HalfEdgeHandle>("visited through", HalfEdgeHandle(-1));

            for(auto out_he: mesh.outgoing_halfedges(v)){
                auto visited_he = visited[mesh.to_vertex_handle(out_he)];
                if(visited_he.idx() == -1){
                    visited[mesh.to_vertex_handle(out_he)] = out_he;
                }else{
                    std::cout<<" FOUND DOUBLE EDGE : "<<std::endl;
                    std::cout<<"    "<<visited_he<<": "<<mesh.halfedge(visited_he)<<std::endl;
                    std::cout<<"    "<<out_he<<": "<<mesh.halfedge(out_he)<<std::endl;
                    return false;
                }
            }
        }

        return true;
    }
}