#include <iostream>
#include <array>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/IO/IO.hh>
#include <OpenVolumeMesh/Unstable/SmartHandles.hh>
#include "tinyusdz.hh"
#include "usda-writer.hh"



namespace OVM = OpenVolumeMesh;

using Mesh = OVM::GeometricPolyhedralMeshV3d;
using VH = OVM::VH;
using CH = OVM::CH;
using FH = OVM::FH;
using EH = OVM::EH;
using HEH = OVM::HEH;
using HFH = OVM::HFH;
using Vec3d = OVM::Vec3d;

void ovm2usd(Mesh const& mesh, tinyusdz::Stage &stage)
{
  std::vector<tinyusdz::value::point3f> vertex_pos;

  std::vector<std::array<double, 3>> vertex_bary;
  //std::vector<std::vector<size_t>> face_indices;
  std::vector<int32_t> face_flat_indices;
  std::vector<int32_t> face_valences;
  vertex_pos.reserve(mesh.n_vertices()*5); // rough guess
  vertex_bary.reserve(mesh.n_vertices()*5); // rough guess
                                            //face_indices.reserve(mesh.n_halffaces());
  face_valences.reserve(mesh.n_halffaces());

  std::vector<int> face_feature_id;
  face_feature_id.reserve(mesh.n_halffaces());

  auto maybe_feature_vert = mesh.get_vertex_property<int>("AlgoHex::FeatureVertices");
  auto maybe_feature_edge = mesh.get_edge_property<int>("AlgoHex::FeatureEdges");
  auto maybe_feature_face = mesh.get_face_property<int>("AlgoHex::FeatureFaces");
  size_t n_faces = 0;
  size_t n_corners = 0;

  auto vh2vi = std::unordered_map<VH, size_t>{};
  for (const auto _ch: mesh.cells()) {
    auto ch = make_smart(_ch, mesh);
    const auto bary = mesh.barycenter(ch);
    vh2vi.clear();
    for (const auto vh: ch.vertices()) {
      const auto &p = mesh.vertex(vh);
      vh2vi[vh] = vertex_pos.size();
      vertex_pos.push_back({(float)p[0], (float)p[1], (float)p[2]});
      vertex_bary.push_back({bary[0], bary[1], bary[2]});
    }
    auto cv = mesh.cell_vertices(_ch);
    for (const auto hfh: ch.halffaces()) {
      size_t n_verts = 0;
      if (maybe_feature_face) {
        face_feature_id.push_back(maybe_feature_face->at(hfh.face_handle()));
      }
      for (auto heh: mesh.halfface_halfedges(hfh)) {
        auto vi = vh2vi[mesh.to_vertex_handle(heh)];
        face_flat_indices.push_back(vi);
        ++n_verts;
      }
      face_valences.push_back(n_verts);
    }
  }

  tinyusdz::Xform xform;
  xform.name = "root";
#if 0
  tinyusdz::XformOp op;
  op.op_type = tinyusdz::XformOp::OpType::Translate;
  tinyusdz::value::double3 translate;
  translate[0] = 1.0;
  translate[1] = 2.0;
  translate[2] = 3.0;
  op.set_value(translate);
  xform.xformOps.push_back(op);
#endif

  tinyusdz::GeomMesh usd_mesh;
  usd_mesh.name = "cells";
  std::cout << "n vert " << vertex_pos.size() <<std::endl;
  std::cout << "n fval " << face_valences.size() <<std::endl;
  std::cout << "n ff " << face_flat_indices.size() <<std::endl;
  usd_mesh.points.set_value(std::move(vertex_pos));
  usd_mesh.faceVertexCounts.set_value(std::move(face_valences));
  usd_mesh.faceVertexIndices.set_value(std::move(face_flat_indices));

  if(1) {
    tinyusdz::Attribute sharp_face_attr;
    std::vector<bool> sharp_face(usd_mesh.get_faceVertexCounts().size(), true);
    tinyusdz::primvar::PrimVar sharp_face_var;
    sharp_face_var.set_value(sharp_face);
    sharp_face_attr.set_var(std::move(sharp_face_var));
    tinyusdz::Property sharp_face_prop(sharp_face_attr, /* custom*/false);
    usd_mesh.props.emplace("primvars:sharp_face", sharp_face_prop);
  }

  tinyusdz::Prim meshPrim(usd_mesh);
  tinyusdz::Prim xformPrim(xform);

  xformPrim.children().emplace_back(std::move(meshPrim));
//stage.add_root_prim(std::move(xformPrim));
  stage.root_prims().emplace_back(std::move(xformPrim));
  //
  //        tinyusdz::Attribute uvAttr;
  //std::vector<tinyusdz::value::texcoord2f> uvs;
  //mesh.props.emplace("primvars:uv", uvProp);
  //
  //
  //std::vector<std::array<double, 3>> meshVertexCellBary;

  //ply.addVertexColors(meshVertexCellBary);
  // ply.addFaceIndices(face_indices);
  // ply.getElement("face").addProperty("face_feature", face_feature_id);
}

int main(int argc, char** argv) {

    Mesh mesh;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <in.ovm(b)> <out.ply>\n";
        return 1;
    }
    const auto path_in = argv[1];
    const auto path_out = argv[2];

    bool success = OVM::IO::read_file(path_in, mesh, false, false);
    if (!success) {
      std::cerr << "Failed to read input mesh." << std::endl;
      return 1;
    }

    // Create an empty object

    tinyusdz::Stage stage;
    ovm2usd(mesh, stage);

    //stage.ExportToString();
    //std::cout << to_string(stage) << "\n";

    std::string warn;
    std::string err;
    bool ret = tinyusdz::usda::SaveAsUSDA("output.usda", stage, &warn, &err);

    if (warn.size()) {
      std::cout << "WARN: " << warn << "\n";
    }

    if (err.size()) {
      std::cerr << "ERR: " << err << "\n";
    }

    //plyOut.write(path_out, happly::DataFormat::Binary);

    return 0;
}
