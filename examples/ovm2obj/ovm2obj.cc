#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/IO/IO.hh>

namespace OVM = OpenVolumeMesh;

using Mesh = OVM::GeometricPolyhedralMeshV3d;
using VH = OVM::VH;
using CH = OVM::CH;
using FH = OVM::FH;
using EH = OVM::EH;
using HEH = OVM::HEH;
using HFH = OVM::HFH;

class SerializeObj {
public:
    SerializeObj(std::ostream &_os, Mesh const&_mesh)
        : os{_os}
        , mesh{_mesh}
    {
        os << std::setprecision(17);

    }
    void serialize() {
        os << "o cells\n";
        n_vertices = 0;
        for (const auto ch: mesh.cells()) {
            serialize_cell(ch);
#if 0
            std::cout << "DEBUG: breaking after 1st cell!" << std::endl;
            break;
#endif
        }
    }
protected:
    size_t serialize_point(VH vh) {
        const auto &p = mesh.vertex(vh);
        os << "v " << p[0] << " "
                   << p[1] << " "
                   << p[2] << "\n";
        return ++n_vertices; //1-based indexing
    }
    void serialize_cell(CH ch) {
        auto vh2vi = std::unordered_map<VH, int>{};
        for (const auto vh: mesh.cell_vertices(ch)) {
            vh2vi[vh] = serialize_point(vh);
        }
        const auto bary = mesh.barycenter(ch);
        os << "vn " << bary[0] << " " << bary[1] << " " << bary[2] << "\n";
        for (const auto hfh: mesh.cell_halffaces(ch)) {
            os << "f";
            // flip orientation on purpose:
            auto opp_hfh = mesh.opposite_halfface_handle(hfh);
            for (const auto heh: mesh.halfface_halfedges(opp_hfh)) {
                auto vi = vh2vi[mesh.to_vertex_handle(heh)];
                os << " " << vi << "//-1";
            }
            os << "\n";
        }

    }
private:
    std::ostream &os;
    Mesh const &mesh;
    size_t n_vertices = 0;
};

int main(int argc, char** argv) {

    Mesh mesh;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <in.ovm(b)> <out.obj>\n";
        return 1;
    }
    const auto path_in = argv[1];
    const auto path_out = argv[2];

    OVM::IO::read_file(path_in, mesh, false, false);

    std::ofstream of(path_out);
    SerializeObj{of, mesh}.serialize();


    return 0;
}
