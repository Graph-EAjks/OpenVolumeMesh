#include "OpenVolumeMesh/Core/TopologyChecks.hh"
#include "OpenVolumeMesh/Core/detail/TopologicalLinkT_impl.hh"
#include <stack>

namespace OpenVolumeMesh {

    CellHandle cell_exists(const TetrahedralMeshTopologyKernel& mesh,
                           const std::vector<VertexHandle>& cell_vertices) {
        if(cell_vertices.size() == 4){
            for(auto vc_it = mesh.vc_iter(cell_vertices[0]); vc_it.is_valid(); ++vc_it){
                int found_vertices_count(0);
                for(auto cv_it = mesh.cv_iter(*vc_it); cv_it.is_valid(); ++cv_it){
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
                                          const VertexHandle vertex,
                                          const FaceHandle face){
        if (!mesh.is_valid(vertex) || !mesh.is_valid(face)) {
            return false;
        }
        auto hf_vertices = mesh.get_halfface_vertices(mesh.halfface_handle(face, 0));

        return  hf_vertices[0] == vertex ||
                hf_vertices[1] == vertex ||
                hf_vertices[2] == vertex;
    }

    bool cell_contains_vertex(const TetrahedralMeshTopologyKernel& mesh,
                                          const VertexHandle vertex,
                                          const CellHandle cell){
        if (!mesh.is_valid(vertex) || !mesh.is_valid(cell)) {
            return false;
        }
        auto c_vertices = mesh.get_cell_vertices(cell);

        return  c_vertices[0] == vertex ||
                c_vertices[1] == vertex ||
                c_vertices[2] == vertex ||
                c_vertices[3] == vertex;

    }

    std::set<std::set<VertexHandle>> find_non_cell_tets(const TetrahedralMeshTopologyKernel& mesh, bool only_check_faces = true){

//        int counter = 0;

        std::set<std::set<VertexHandle>> non_cell_tets;

        // TODO: Maybe it can be faster by marking some vertices/triangles as visited?
        if (!only_check_faces) {
            auto non_face_tris = OpenVolumeMesh::find_non_face_triangles(mesh);
            // check for all non-face triangles, if there exists a fourth vertex which forms a tet with the triangle
            for (auto triangle : non_face_tris) {
                auto potential_tet_vertex = mesh.create_private_property<int, Entity::Vertex>("potential_tet_vertex");
                std::vector<VertexHandle> triangle_vertices;
                for (auto tri_vertex : triangle) {
                    triangle_vertices.push_back(tri_vertex);
                }

                // iterate over the neighbors of tri_vertex0 and mark them as a potential tet extension
                for (auto vv_it = mesh.vv_iter(triangle_vertices[0]); vv_it.is_valid(); ++vv_it) {
                    ++potential_tet_vertex[*vv_it];
                }
                // iterate over the neighbors of tri_vertex1 and mark them as a potential tet extension
                for (auto vv_it = mesh.vv_iter(triangle_vertices[1]); vv_it.is_valid(); ++vv_it) {
                    ++potential_tet_vertex[*vv_it];
                }
                // iterate over the neighbors of tri_vertex2. If any of them is a neighbor of the first two vertices as well,
                // then we found a non-cell tet
                for (auto vv_it = mesh.vv_iter(triangle_vertices[2]); vv_it.is_valid(); ++vv_it) {
                    if (potential_tet_vertex[*vv_it] == 2){
                        std::set<VertexHandle> non_cell_tet;

                        non_cell_tet.insert(triangle_vertices[0]);
                        non_cell_tet.insert(triangle_vertices[1]);
                        non_cell_tet.insert(triangle_vertices[2]);
                        non_cell_tet.insert(*vv_it);
                        non_cell_tets.insert(non_cell_tet);
                    }
                }
            }
        }

        for(auto v : mesh.vertices()){

//            auto neighbor = mesh.request_vertex_property<bool>("neighbor");
            auto neighbor = mesh.create_private_property<bool, Entity::Vertex>("neighbor");

            for(auto vv_it = mesh.vv_iter(v); vv_it.is_valid(); ++vv_it){
                neighbor[*vv_it] = true;
            }

            for(auto vf_it = mesh.vf_iter(v); vf_it.is_valid(); ++vf_it){

                std::set<VertexHandle> tet_vertices;
                tet_vertices.insert(v);


                //find opposite edge to v
                EdgeHandle opposite_edge(-1);

                for(auto fe_it = mesh.fe_iter(*vf_it); fe_it.is_valid(); ++fe_it){
                    if(mesh.edge(*fe_it).from_vertex() != v &&
                       mesh.edge(*fe_it).to_vertex() != v){
                        tet_vertices.insert(mesh.edge(*fe_it).from_vertex());
                        tet_vertices.insert(mesh.edge(*fe_it).to_vertex());

                        opposite_edge = *fe_it;
                    }
                }

                //iterate through the opposite edge's faces to try to find another neighbor vertex
                for(auto ef_it = mesh.ef_iter(opposite_edge); ef_it.is_valid(); ++ef_it){
                    std::set<VertexHandle> partial_cell = tet_vertices;

                    //find the opposite vertex
                    for(auto fv_it = mesh.fv_iter(*ef_it); fv_it.is_valid(); ++fv_it){
                        if(mesh.edge(opposite_edge).from_vertex() != *fv_it &&
                           mesh.edge(opposite_edge).to_vertex() != *fv_it &&
                           neighbor[*fv_it]){

                            partial_cell.insert(*fv_it);

                            bool found_cell(false);

                            //check if it's a cell or not
                            for(auto vc_it = mesh.vc_iter(v); vc_it.is_valid(); ++vc_it){
                                int found_vertices(0);
                                for(auto cv_it = mesh.cv_iter(*vc_it); cv_it.is_valid(); ++cv_it){
                                    if(partial_cell.find(*cv_it) != partial_cell.end()){
                                        found_vertices++;
                                    }
//                                    ++counter;
                                }
                                if(found_vertices == 4){
                                    found_cell = true;
                                }
                            }
                            bool faces_exist = true;
                            std::vector<VertexHandle> cell_vertices;
                            for (auto v2 : partial_cell) {
                                cell_vertices.push_back(v2);
                            }
                            std::vector<std::vector<VertexHandle>> triangles = {{cell_vertices[0], cell_vertices[1], cell_vertices[2]},
                                                                                {cell_vertices[0], cell_vertices[1], cell_vertices[3]},
                                                                                {cell_vertices[0], cell_vertices[2], cell_vertices[3]},
                                                                                {cell_vertices[1], cell_vertices[2], cell_vertices[3]}};
                            for (auto triangle : triangles) {
                                faces_exist &= mesh.halfface(triangle).is_valid();
                            }

                            if(!found_cell && (!only_check_faces || faces_exist)){

                                if(partial_cell.size() != 4){
                                    throw Core::detail::invalid_topology_error("tet does not contain 4 vertices");
                                    return {};//TODO: is this dead code? I think so, but not sure
                                }
                                non_cell_tets.insert(partial_cell);
                            }
                        }
                    }
                }
            }
        }
//        std::cout << "Version 1: " << counter << std::endl;
        return non_cell_tets;
    }

    std::set<std::set<VertexHandle>> find_non_cell_tets_2(const TetrahedralMeshTopologyKernel& mesh, bool only_check_faces) {

//        int counter = 0;

        std::set<std::set<VertexHandle>> ret;

        auto visited = mesh.create_private_property<bool, Entity::Vertex>("visited", false);

        for(auto v: mesh.vertices()){
            visited[v] = true;
            for (auto vv_it = mesh.vv_iter(v); vv_it.is_valid(); ++vv_it) {
                if (visited[*vv_it]) continue;
                for (auto vv_it_2 = mesh.vv_iter(v); vv_it_2.is_valid(); ++vv_it_2) {
                    if (*vv_it_2 == *vv_it) continue;
                    if (visited[*vv_it_2]) continue;
                    // first, we check, if vv_it_2 is a neighbour of vv_it.
                    bool found = false;
                    for (auto vv_it_3 = mesh.vv_iter(*vv_it_2); vv_it_3.is_valid(); ++vv_it_3) {
                        if (*vv_it_3 == *vv_it) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        break; // vv_it_2 is not a neighbour of vv_it
                    }
                    // we have found a triangle, now let's try to find another vertex to make it a tet
                    for (auto vv_it_3 = mesh.vv_iter(v); vv_it_3.is_valid(); ++vv_it_3) {
                        if (*vv_it_3 == *vv_it || *vv_it_3 == *vv_it_2) continue;
                        // first, we check, if vv_it_3 is a neighbour of vv_it.
                        found = false;
                        for (auto vv_it_4 = mesh.vv_iter(*vv_it_3); vv_it_4.is_valid(); ++vv_it_4) {
                            if (*vv_it_4 == *vv_it) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            break; //v_it_3 is not a neighbour of vv_it
                        }
                        // second, we check, if vv_it_3 is a neighbour of vv_it_2
                        found = false;
                        for (auto vv_it_4 = mesh.vv_iter(*vv_it_3); vv_it_4.is_valid(); ++vv_it_4) {
                            if (*vv_it_4 == *vv_it_2) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            break; //vv_it_3 is not a neighbour of vv_it_2
                        }
                        // now, we know, that v, vv_it, v_it_2 and vv_it_3 form a tet, now we need to check if the cell exists
                        bool found_cell = false;
                        for (auto vc_it = mesh.vc_iter(v); vc_it.is_valid(); ++vc_it) {
                            int vertex_count = 0;
                            for (auto cv_it = mesh.cv_iter(*vc_it); cv_it.is_valid(); ++cv_it) {
                                if (*cv_it == *vv_it ||
                                    *cv_it == *vv_it_2 ||
                                    *cv_it == *vv_it_3) {
                                    ++vertex_count;
                                }
//                                ++counter;
                            }
                            if (vertex_count == 3) {
                                found_cell = true;
                                break;
                            }
                        }
                        bool faces_exist = true;
                        if (only_check_faces) {
                            std::vector<std::vector<VertexHandle>> triangles = {{v,*vv_it, *vv_it_2},
                                                                                {v,*vv_it, *vv_it_3},
                                                                                {v,*vv_it_2, *vv_it_3},
                                                                                {*vv_it,*vv_it_2, *vv_it_3}};
                            for (auto triangle : triangles) {
                                faces_exist &= mesh.halfface(triangle).is_valid();
                            }
                        }
                        if (!found_cell && faces_exist) {
                            std::set<VertexHandle> non_cell_tet;
                            non_cell_tet.insert(v);
                            non_cell_tet.insert(*vv_it);
                            non_cell_tet.insert(*vv_it_2);
                            non_cell_tet.insert(*vv_it_3);
                            ret.insert(non_cell_tet);
                        }
                    }
                }
            }
        }
//        std::cout << "Version 2: " << counter << std::endl;
        return ret;
    }

    bool link_condition(const TetrahedralMeshTopologyKernel& mesh,
                        const EdgeHandle edge){

        auto from_vertex = mesh.edge(edge).from_vertex();
        auto to_vertex = mesh.edge(edge).to_vertex();

        auto vertex_link_intersection = link(mesh, from_vertex).intersection(link(mesh, to_vertex));

        if(vertex_link_intersection.faces().size()){
            return false;
        }

        return vertex_link_intersection == link(mesh, edge);
    }

    bool single_connected_component(const TetrahedralMeshTopologyKernel&  mesh){

        if (mesh.n_logical_vertices() == 0) {// mesh contains no vertices, so it is empty
            return false;
        }

        auto visited_prop = mesh.create_private_property<bool, Entity::Vertex>("visited", false);

//        std::set<VertexHandle> to_visit;
        std::stack<VertexHandle> to_visit;
//        to_visit.insert(*mesh.vertices().first);
        to_visit.push(*mesh.vertices().first);

        while(!to_visit.empty()){
//            auto next = *to_visit.begin();
            auto next = to_visit.top();
//            to_visit.erase(to_visit.begin());
            to_visit.pop();

            visited_prop[next] = true;

            for(auto out_he_it = mesh.outgoing_halfedges(next).first; out_he_it.is_valid(); ++out_he_it){
                auto to_vertex = mesh.to_vertex_handle(*out_he_it);
                if(!visited_prop[to_vertex]){
//                    to_visit.insert(to_vertex);
                    to_visit.push(to_vertex);
                }
            }
        }

        bool all_visited(true);

        //TODO: count visited to be faster?
        for(auto v: mesh.vertices()){
            if(!visited_prop[v]){
                all_visited = false;
            }
        }

        return all_visited;
    }

    // largely copied from OpenFlipper-Free::Plugin-InfoVolumeMeshObject::VolumeMeshAnalysis
    size_t count_connected_components(const TetrahedralMeshTopologyKernel& mesh) {
        auto visited = mesh.create_private_property<bool, Entity::Vertex>("visited", false);
        size_t n_connected_components = 0;

        size_t n_vertices_found = 0; // for early abort
        size_t n_vertices = mesh.n_logical_vertices();

        for (const auto &root_vh: mesh.vertices())
        {
            if(n_vertices_found == n_vertices)
                break;
            if (visited[root_vh])
                continue;
            ++n_connected_components;
            visited[root_vh] = true;
            ++n_vertices_found;

            std::stack<VH> stack;
            stack.push(root_vh);
            while(!stack.empty())
            {
                const auto vh = stack.top(); stack.pop();
                for (const auto &nb: mesh.vertex_vertices(vh))
                {
                    if (!visited[nb]) {
                        visited[nb] = true;
                        stack.push(nb);
                        ++n_vertices_found;
                    }
                }
            }
        }
        return n_connected_components;
    }

    bool contains_void(const TetrahedralMeshTopologyKernel&  mesh){

        if (!mesh.vertices_begin().is_valid()) {
            // mesh is empty
            return false;
        }

//        std::vector<FaceHandle> to_visit;
        std::stack<FaceHandle> to_visit;

        for(auto f: mesh.faces()){
            if(mesh.is_boundary(f)){
//                to_visit.push_back(f);
                to_visit.push(f);
                break;
            }
        }

        if(!to_visit.size()){
            throw Core::detail::invalid_topology_error("There is no boundary face, which is impossible");
        }

        auto visited = mesh.create_private_property<bool, Entity::Face>("visited", false);

        int iteration_count(0);

        while(to_visit.size()){
            //std::cout<<" stack size at iteration "<<iteration_count<<": "<<to_visit.size()<<std::endl;

//            auto current_face = to_visit.back();
            auto current_face = to_visit.top();
//            to_visit.pop_back();
            to_visit.pop();

            if(!visited[current_face]){

                visited[current_face] = true;

                for(auto fe_it = mesh.fe_iter(current_face); fe_it.is_valid(); ++fe_it){

                    for(auto ef_it = mesh.ef_iter(*fe_it); ef_it.is_valid(); ++ef_it){
                        if(mesh.is_boundary(*ef_it) && !visited[*ef_it]){
//                            to_visit.push_back(*ef_it);
                            to_visit.push(*ef_it);
                        }
                    }
                }
            }

            iteration_count++;
        }


        for(auto f: mesh.faces()){
            if(mesh.is_boundary(f) && !visited[f]){
//                std::cout<<" found a face not visited"<<std::endl;
                return true;
            }
        }

        return false;
    }

    bool manifold_vertex(const TetrahedralMeshTopologyKernel& mesh,
                         const VertexHandle vertex){


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
        for(auto vf_it = mesh.vf_iter(vertex); vf_it.is_valid(); ++vf_it){
            if(mesh.is_boundary(*vf_it)){
                initial_boundary_face = *vf_it;
                break;
            }
        }

        if(initial_boundary_face.idx() == -1){
//            std::cerr<<" couldn't find boundary face adjacent to boundary vertex. returning false"<<std::endl;
            return false;
        }

        //std::cout<<" ----------------------------"<<std::endl;
        //std::cout<<" checking vertex "<<vertex<<std::endl;

        //std::cout<<" -- initial boundary face : "<<initial_boundary_face<<": "<<mesh.face(initial_boundary_face)<<std::endl;

        auto visited_edge = mesh.create_private_property<bool, Entity::Edge>("visited", false);
        auto visited_face = mesh.create_private_property<bool, Entity::Face>("visited", false);

        EdgeHandle next_edge(-1);

        //visit all vertices of the initial face
        for(auto fe_it = mesh.fe_iter(initial_boundary_face); fe_it.is_valid(); ++fe_it){
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
            for(auto ef_it = mesh.ef_iter(current_edge); ef_it.is_valid(); ++ef_it){
                //std::cout<<" -- checking face "<<*ef_it<<": "<<mesh.face(*ef_it)<<std::endl;
                if(mesh.is_boundary(*ef_it) && ! visited_face[*ef_it]){
                    //found a neighboring face to both current_edge and vertex
                    //so visit the next edge and use it as next current_edge

                    //std::cout<<" ---- neighbor face: "<<*ef_it<<": "<<mesh.face(*ef_it)<<std::endl;

                    for(auto fe_it = mesh.fe_iter(*ef_it); fe_it.is_valid(); ++fe_it){
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
        for(auto vf_it = mesh.vf_iter(vertex); vf_it.is_valid(); ++vf_it){
            if(mesh.is_boundary(*vf_it)){
                all_visited &= visited_face[*vf_it];
            }
        }
//        if(!all_visited){
//            std::cout<<" vertex "<<vertex<<" is non-manifold"<<std::endl;
//        }

        return all_visited;
    }

    std::optional<std::pair<HEH, HEH>> contains_double_edges(const TetrahedralMeshTopologyKernel& mesh) {
        for(auto v: mesh.vertices()){
            // TODO: replace with smart tagger -> complete unittests for smart tagger first?
            auto visited = mesh.create_private_property<HalfEdgeHandle, Entity::Vertex>("visited through", HalfEdgeHandle(-1));

            for(auto out_he: mesh.outgoing_halfedges(v)){
                auto visited_he = visited[mesh.to_vertex_handle(out_he)];
                if(visited_he.idx() == -1){
                    visited[mesh.to_vertex_handle(out_he)] = out_he;
                }else{
                    std::pair<HEH, HEH> ret;
                    ret.first = visited_he;
                    ret.second = out_he;
                    return ret;
                }
            }
        }
        return std::nullopt;
    }

    /*
    bool no_double_edges(const TetrahedralMeshTopologyKernel& mesh){

        //TODO: return std::optional<std::pair<HEH, HEH>> which contains exactly one example of a double edge or nothing instead of bool

        for(auto v: mesh.vertices()){
            auto visited = mesh.create_private_property<HalfEdgeHandle, Entity::Vertex>("visited through", HalfEdgeHandle(-1));

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
     */

    std::set<std::set<HEH>> find_multi_edges(const TetrahedralMeshTopologyKernel& mesh) {

        std::set<std::set<HEH>> multi_edges;
        for(auto v: mesh.vertices()){
            std::set<HEH> empty_heh_set;
            auto multi_edge_list = mesh.create_private_property<std::set<HEH>, Entity::Vertex>("multi_edge_list", empty_heh_set);

            for(auto out_he: mesh.outgoing_halfedges(v)){
                multi_edge_list[mesh.to_vertex_handle(out_he)].insert(out_he);
            }
            for (auto vv_it = mesh.vv_iter(v); vv_it.is_valid(); ++vv_it) {
                if (multi_edge_list[*vv_it].size() > 1) {
                    multi_edges.insert(multi_edge_list[*vv_it]);
                }
            }
        }
        return multi_edges;
    }

    /*
    void print_mesh_topology(const TetrahedralMeshTopologyKernel& mesh) {
        std::cout << "Printing mesh topology: " << std::endl;
        std::cout << std::endl;
        std::cout << "Number ob vertices: " << mesh.n_vertices() << std::endl;
        std::cout << std::endl;
        std::cout << "Number of edges: " << mesh.n_edges() << std::endl;
        int counter = 0;
        for (auto e_it = mesh.edges_begin(); e_it.is_valid(); ++e_it) {
            auto heh = mesh.halfedge_handle(*e_it, 0);
            std::cout << counter << ": " << mesh.from_vertex_handle(heh).idx() << "-" << mesh.to_vertex_handle(heh).idx() << std::endl;
            ++counter;
        }
        std::cout << std::endl;
        std::cout << "Number of faces: " << mesh.n_faces() << std::endl;
        counter = 0;
        for (auto f_it = mesh.faces_begin(); f_it.is_valid(); ++f_it) {
            std::cout << counter << ": ";
            for (auto f_he_it = mesh.fhe_iter(*f_it); f_he_it.is_valid(); ++f_he_it) {
                std::cout << (*f_he_it).idx() << "-";
            }
            std::cout << '\b' << " " << std::endl;
            ++counter;
        }
        std::cout << std::endl;
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
        counter = 0;
        for (auto c_it = mesh.cells_begin(); c_it.is_valid(); ++c_it) {
            std::cout << counter << ": ";
            for (auto c_hf_it = mesh.chf_iter(*c_it); c_hf_it.is_valid(); ++c_hf_it) {
                std::cout << (*c_hf_it).idx() << "-";
            }
            std::cout << '\b' << " " << std::endl;
            ++counter;
        }
    }
     */

    std::ostream &operator<<(std::ostream &os, const TetrahedralMeshTopologyKernel &mesh) {
        os  << "Printing mesh topology: " << std::endl;
        std::cout << std::endl;
        std::cout << "Number ob vertices: " << mesh.n_vertices() << std::endl;
        std::cout << std::endl;
        std::cout << "Number of edges: " << mesh.n_edges() << std::endl;
        int counter = 0;
        for (auto e_it = mesh.edges_begin(); e_it.is_valid(); ++e_it) {
            auto heh = mesh.halfedge_handle(*e_it, 0);
            std::cout << counter << ": " << mesh.from_vertex_handle(heh).idx() << "-" << mesh.to_vertex_handle(heh).idx() << std::endl;
            ++counter;
        }
        std::cout << std::endl;
        std::cout << "Number of faces: " << mesh.n_faces() << std::endl;
        counter = 0;
        for (auto f_it = mesh.faces_begin(); f_it.is_valid(); ++f_it) {
            std::cout << counter << ": ";
            for (auto f_he_it = mesh.fhe_iter(*f_it); f_he_it.is_valid(); ++f_he_it) {
                std::cout << (*f_he_it).idx() << "-";
            }
            std::cout << '\b' << " " << std::endl;
            ++counter;
        }
        std::cout << std::endl;
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
        counter = 0;
        for (auto c_it = mesh.cells_begin(); c_it.is_valid(); ++c_it) {
            std::cout << counter << ": ";
            for (auto c_hf_it = mesh.chf_iter(*c_it); c_hf_it.is_valid(); ++c_hf_it) {
                std::cout << (*c_hf_it).idx() << "-";
            }
            std::cout << '\b' << " " << std::endl;
            ++counter;
        }
        return os;
    }

    std::set<std::set<VertexHandle>> find_non_face_triangles(const TetrahedralMeshTopologyKernel& mesh) {
        std::set<std::set<VertexHandle>> found_faces;
        auto visited = mesh.create_private_property<bool, Entity::Vertex>("visited", false);
        for (auto v_it = mesh.vertices_begin(); v_it.is_valid(); ++v_it) {
            auto vertex = *v_it;
            auto found_faces_by_vertex = find_non_face_triangles_around_vertex(mesh, vertex);
            found_faces.insert(found_faces_by_vertex.begin(), found_faces_by_vertex.end());
            visited[*v_it] = true;
        }
        return found_faces;
    }

    std::set<std::set<VertexHandle>> find_non_face_triangles_around_vertex(const TetrahedralMeshTopologyKernel& mesh, const VertexHandle vertex) {
        std::set<std::set<VertexHandle>> res;

        // call visited property from find_non_face_triangles method to prevent unnecessary checks
        auto visited = mesh.create_private_property<bool, Entity::Vertex>("visited", false);

        // mark the neighbors
        auto neighbor_prop = mesh.create_private_property<bool, Entity::Vertex>("neighbor_prop", false);
        std::set<VertexHandle> neighbors;
        for (auto vv_it = mesh.vv_iter(vertex); vv_it.is_valid(); ++vv_it) {
            if (!visited[*vv_it]) {
                neighbor_prop[*vv_it] = true;
                neighbors.insert(*vv_it);
            }
        }

        for (auto neighbor : neighbors) {
            for (auto vv_it = mesh.vv_iter(neighbor); vv_it.is_valid(); ++vv_it) {
                if (neighbor_prop[*vv_it]) {
                    std::set<VertexHandle> found_triangle = {vertex, neighbor, *vv_it};
                    std::vector<VertexHandle> triangle_vector(found_triangle.size());
                    std::copy(found_triangle.begin(), found_triangle.end(), triangle_vector.begin());
                    if (!mesh.halfface(triangle_vector).is_valid()) {
                        res.insert(found_triangle);
                    }
                }
            }
        }
        return res;
    }
}