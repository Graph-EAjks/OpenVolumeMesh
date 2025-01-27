#include <iostream>
#include <array>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/IO/IO.hh>
#include <OpenVolumeMesh/Unstable/SmartHandles.hh>
#include <happly.h>

namespace OVM = OpenVolumeMesh;

using Mesh = OVM::GeometricPolyhedralMeshV3d;
using VH = OVM::VH;
using CH = OVM::CH;
using FH = OVM::FH;
using EH = OVM::EH;
using HEH = OVM::HEH;
using HFH = OVM::HFH;
using Vec3d = OVM::Vec3d;

#if 0
class PLYFile {
    void serialize() {


        os << "ply\nformat ascii 1.0\n";
        auto n_corners = mesh.n_cells() * 8; // XXX TODO for hex meshes
        auto n_faces = mesh.n_cells() * 6; // XXX TODO for hex meshes
        os << "element vertex " << n_corners << "\n";
        os << "property float x\n"
           << "property float y\n"
           << "property float z\n"
           << "property float bary_x\n"
           << "property float bary_y\n"
           << "property float bary_z\n"
           << "property float ch\n"
           ;
        os << "element face " << n_faces << "\n"
           << "property list uchar int vertex_indices\n"
           << "property float feat\n";
        os << "element edge " << n_edges << "\n"
           << "property list uchar int vertex_indices\n"
           << "property float feat\n";
        os << "end_header\n";
        faces.reserve(n_faces);
        n_vertices = 0;
        for (const auto ch: mesh.cells()) {
            serialize_cell(ch);
#if 0
            std::cout << "DEBUG: breaking after 1st cell!" << std::endl;
            break;
#endif
        }
        for (const auto &[verts, f_id]: faces) {
            os << verts.size();
            for (const auto v: verts) {
                os << " " << v;
            }
            os << " " << f_id;
            os << "\n";
        }
    }
protected:
    struct Point {
        Vec3d pos;
        Vec3d bary;
        CH ch;
    };
    size_t serialize_point(VH vh, OVM::Vec3d const& bary, CH ch) {
        const auto &p = mesh.vertex(vh);
        os << p[0] << " "
           << p[1] << " "
           << p[2] << " "
           << bary[0] << " "
           << bary[1] << " "
           << bary[2] << " "
           << ch.idx() << "\n";
        return n_vertices++;
    }
    void serialize_cell(CH ch) {
        const auto bary = mesh.barycenter(ch);
        auto vh2vi = std::unordered_map<VH, int>{};
        for (const auto vh: mesh.cell_vertices(ch)) {
            vh2vi[vh] = serialize_point(vh, bary, ch);
        }
        std::vector<int> face;
        for (const auto hfh: mesh.cell_halffaces(ch)) {
            auto f_id = feature_face[hfh.face_handle()];
            face.clear();
            // flip orientation on purpose:
            auto opp_hfh = mesh.opposite_halfface_handle(hfh);
            for (const auto heh: mesh.halfface_halfedges(opp_hfh)) {
                auto vi = vh2vi[mesh.to_vertex_handle(heh)];
                face.push_back(vi);
            }
            faces.emplace_back(std::move(face), f_id);
        }

    }

private:
    std::ofstream os;
    Mesh const &mesh;
    OVM::PropertyPtr<int,OVM::Entity::Face> feature_face = *mesh.get_face_property<int>("AlgoHex::FeatureFaces");
    std::vector<std::pair<std::vector<int>, int>> faces; // XXX todo 4 for hex
    size_t n_vertices = 0;
};
class SerializeCellsPLY {
public:
    SerializeCellsPLY(std::filesystem::path const &path, Mesh const&_mesh)
        : os{path.c_str()}
        , mesh{_mesh}
    {
        os << std::setprecision(17);
    }


class SerializeFeaturesOBJ {
public:
    SerializeFeaturesOBJ(std::filesystem::path const &path, Mesh const&_mesh)
        : os{path.c_str()}
        , mesh{_mesh}
    {
        os << std::setprecision(17);

    }
    void serialize() {
        if (auto maybe_feature_vert = mesh.get_vertex_property<int>("AlgoHex::FeatureVertices")) {
            auto feature_vert = *maybe_feature_vert;
            new_object("feature_verts");
            for (const auto vh: mesh.vertices()) {
                if (feature_vert[vh]) {
                    serialize_point(vh);
                }
            }
        }
        if (auto maybe_feature_edge = mesh.get_edge_property<int>("AlgoHex::FeatureEdges")) {
            new_object("feature_edges");
            auto feature_edge = *maybe_feature_edge;
            for (const auto eh: mesh.edges()) {
                auto f_id = feature_edge[eh];
                if (!f_id) continue;
                auto v_from = mesh.to_vertex_handle(eh.halfedge_handle(0));
                auto v_to = mesh.to_vertex_handle(eh.halfedge_handle(1));
                auto vi0 = serialize_point(v_from);
                auto vi1 = serialize_point(v_to);
                os << "vt " << f_id << " " << 0 << "\n";
                //os << "l " << vi0 << "/-1 " << vi1 << "/-1\n";
                os << "l " << vi0 << "/-1 " << vi1 << "/-1\n";
            }
        }
    }
protected:
    void new_object(std::string const &name) {
        std::cout << name << std::endl;
        os << "o " << name << "\n";
        //n_vertices = 0;
        //vi.fill(0);
    }
    size_t serialize_point(VH vh) {
        if (vi[vh] > 0) return vi[vh];
        const auto &p = mesh.vertex(vh);
        os << "v " << p[0] << " "
                   << p[1] << " "
                   << p[2] << "\n";
        auto vidx = ++n_vertices; //1-based indexing
        vi[vh] = vidx;
        return vidx;
    }
private:
    std::ofstream os;
    Mesh const &mesh;
    OVM::PropertyPtr<int,OVM::Entity::Vertex> vi = mesh.create_private_property<int, OVM::Entity::Vertex>();
    size_t n_vertices = 0;
};
#endif

void ovm2ply(Mesh const& mesh, happly::PLYData &ply)
{
    std::vector<std::array<double, 3>> vertex_pos;
    std::vector<std::array<double, 3>> vertex_bary;
    std::vector<std::vector<size_t>> face_indices;
    vertex_pos.reserve(mesh.n_vertices()*5); // rough guess
    vertex_bary.reserve(mesh.n_vertices()*5); // rough guess
    face_indices.reserve(mesh.n_halffaces());
    std::vector<int> face_feature_id;
    face_feature_id.reserve(mesh.n_halffaces());

    auto maybe_feature_vert = mesh.get_vertex_property<int>("AlgoHex::FeatureVertices");
    auto maybe_feature_edge = mesh.get_edge_property<int>("AlgoHex::FeatureEdges");
    auto maybe_feature_face = mesh.get_face_property<int>("AlgoHex::FeatureFaces");
    size_t n_faces = 0;
    size_t n_verts = 0;
    size_t n_corners = 0;

    auto vh2vi = std::unordered_map<VH, size_t>{};
    for (const auto _ch: mesh.cells()) {
        auto ch = make_smart(_ch, mesh);
        const auto bary = mesh.barycenter(ch);
        vh2vi.clear();
        for (const auto vh: ch.vertices()) {
            const auto &p = mesh.vertex(vh);
            vh2vi[vh] = vertex_pos.size();
            vertex_pos.push_back({p[0], p[1], p[2]});
            vertex_bary.push_back({bary[0], bary[1], bary[2]});
        }
        auto cv = mesh.cell_vertices(_ch);
        for (const auto hfh: ch.halffaces()) {
            face_indices.emplace_back();
            auto &f = face_indices.back();
            f.reserve(mesh.valence(hfh.face_handle()));
            if (maybe_feature_face) {
                face_feature_id.push_back(maybe_feature_face->at(hfh.face_handle()));
            }
            for (auto heh: mesh.halfface_halfedges(hfh)) {
                auto vi = vh2vi[mesh.to_vertex_handle(heh)];
                f.push_back(vi);
            }
        }
    }
    //std::vector<std::array<double, 3>> meshVertexCellBary;

    ply.addVertexPositions(vertex_pos);
    //ply.addVertexColors(meshVertexCellBary);
    ply.addFaceIndices(face_indices);
    ply.getElement("face").addProperty("face_feature", face_feature_id);
}

int main(int argc, char** argv) {

    Mesh mesh;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <in.ovm(b)> <out.ply>\n";
        return 1;
    }
    const auto path_in = argv[1];
    const auto path_out = argv[2];

    OVM::IO::read_file(path_in, mesh, false, false);

    // Create an empty object
    happly::PLYData plyOut;

    ovm2ply(mesh, plyOut);

    plyOut.write(path_out, happly::DataFormat::ASCII);
    //plyOut.write(path_out, happly::DataFormat::Binary);

    return 0;
}
