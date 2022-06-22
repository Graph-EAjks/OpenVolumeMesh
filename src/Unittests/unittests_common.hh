#pragma once

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/FileManager/TypeNames.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

#ifdef __clang__
#  pragma GCC diagnostic ignored "-Weverything"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wglobal-constructors"
#  pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#  pragma GCC diagnostic ignored "-Wmissing-noreturn"
#endif

#include <gtest/gtest.h>

#define EXPECT_HANDLE_EQ(a, b)  EXPECT_EQ((a).idx(), (b).idx())
#define EXPECT_HANDLE_NE(a, b)  EXPECT_NE((a).idx(), (b).idx())


/*
 * Simple test setting for polyhedral meshes
 */

typedef OpenVolumeMesh::GeometricPolyhedralMeshV3d PolyhedralMesh;

class PolyhedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generatePolyhedralMesh(PolyhedralMesh& _mesh);

  // This member will be accessible in all tests
  PolyhedralMesh mesh_;
};

/*
 * Simple test setting for hexahedral meshes
 */

typedef OpenVolumeMesh::GeometricHexahedralMeshV3d HexahedralMesh;

class HexahedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generateHexahedralMesh(HexahedralMesh& _mesh);

  // This member will be accessible in all tests
  HexahedralMesh mesh_;
};


/*
 * Simple test setting for tetrahedral meshes
 */

typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d TetrahedralMesh;

class TetrahedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic tetrahedral mesh consisting of one tet
  void generateTetrahedralMesh(TetrahedralMesh& _mesh);
  // Generate a basic tetrahedral mesh consisting of two tets sharing exactly one face
  void generateTetrahedralMesh_2(TetrahedralMesh& _mesh);
  // generate a basic tetrahedral mesh consisting of one tet. However, the tet cell does not exist, while all
  // vertices, edges and faces do exist
  void generateTetWithoutCell(TetrahedralMesh &mesh);
  // Generate a basic tetrahedral mesh consisting of one triangle, which is not a face
  void generateTriWithoutFace(TetrahedralMesh& mesh);
  // Generate a tetrahedral mesh consisting of two tets without cells or faces
  void generateTet_withoutCellsAndFaces(TetrahedralMesh& mesh);

  //-------------------------------//
  // Non-manifold test tets. Note: They have no geometry, only topology
  //-------------------------------//
  // Generate a non-manifold tetrahedral mesh, consisting of two tets, sharing exactly one vertex
  void generateNonManifoldTet_2T1V(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of two tets, sharing exactly one edge
  void generateNonManifoldTet_2T1E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of two fans of cells sharing exactly one edge
  void generateNonManifoldTet_6T1E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of one tet and one face sharing exactly one edge
  void generateNonManifoldTet_1T1F1E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of two faces, sharing exactly one edge
  void generateNonManifoldTet_2F1E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of three tets sharing one vertex among all of them
  // and each pair sharing one edge. The three shared edges form a triangle. The three tets form a fourth tet,
  // which is no cell
  void generateNonManifoldTet_3T1V3E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of a single vertex
  void generateNonManifoldTet_1V(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of a single edge
  void generateNonManifoldTet_1E(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of a single face
  void generateNonManifoldTet_1F(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of one tet and one face sharing exactly one vertex
  void generateNonManifoldTet_1T1F1V(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of two faces, sharing exactly one vertex
  void generateNonManifoldTet_2F1V(TetrahedralMesh& mesh);
  // Generate a non-manifold tetrahedral mesh, consisting of one face and one edge, sharing exactly one vertex
  void generateNonManifoldTet_1F1E1V(TetrahedralMesh& mesh);

  // This member will be accessible in all tests
  TetrahedralMesh mesh_;
};


// Printer class (for STL compliance test)
class Print {
public:
  explicit Print(bool _mute = false) : mute_(_mute) {}
  void mute(bool _mute) { mute_ = _mute; }
  bool mute() const { return mute_; }
  template<typename Handle>
  void operator()(const Handle &_h) const {
    if(!mute_) std::cerr << "Handle: " << _h << std::endl;
  }
private:
  bool mute_;
};

