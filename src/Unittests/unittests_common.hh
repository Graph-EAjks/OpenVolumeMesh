#pragma once

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/FileManager/TypeNames.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Core/detail/TopologicalFaceSet.hh>

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
  static void generate_tetrahedral_mesh(TetrahedralMesh& _mesh);
  // Generate a basic tetrahedral mesh consisting of two tets sharing exactly one face
  void generate_tetrahedral_mesh_2(TetrahedralMesh& _mesh);
  // generate a basic tetrahedral _mesh consisting of one tet. However, the tet cell does not exist, while all
  // vertices, edges and faces do exist
  void generate_tet_without_cell(TetrahedralMesh& _mesh);
  // Generate a basic tetrahedral _mesh consisting of one triangle, which is not a face
  void generate_tri_without_face(TetrahedralMesh& _mesh);
  // Generate a tetrahedral _mesh consisting of two tets without cells or faces
  void generate_tet_without_cells_and_faces(TetrahedralMesh& _mesh);
  // Generate a tritet: A tetrahedral _mesh basically consisting of two tets, sharing exactly one face, where the two
  // vertices not belonging to that face are connected by an edge, leading to three tets, pairwise sharing one face
  // (contains only vertices and edges, no faces or cells
  void generate_tri_tet(TetrahedralMesh& _mesh);
  // Generate the same as before, but with faces
  void generate_tri_tet_with_faces(TetrahedralMesh& _mesh);
  // Generate a nested tetrahedral _mesh consisting of two similar tets, where the inner and outer tets are connected by
  // edges connecting one vertex of the outer tet with one vertex of the inner tet. This leads to two tets and three
  // polyhedra with two triangles and three quadrilateral faces. Vertices, edges and faces exists, cells do not.
  void generate_nested_tets(TetrahedralMesh& _mesh);
  // Generate a tet with with only three faces
  void generate_tet_3F(TetrahedralMesh& _mesh);
  // Genrete two tets, which are not connected, so there are two connected components in the mesh
  void generate_tets_two_connected_components(TetrahedralMesh& _mesh);

  //-------------------------------//
  // Non-manifold test tets. Note: They have no geometry, only topology
  //-------------------------------//
  // Generate a non-manifold tetrahedral _mesh, consisting of two tets, sharing exactly one vertex
  void generate_non_manifold_tet_2T1V(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of two tets, sharing exactly one edge
  void generate_non_manifold_tet_2T1E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of two fans of cells sharing exactly one edge
  void generate_non_manifold_tet_6T1E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of one tet and one face sharing exactly one edge
  void generate_non_manifold_tet_1T1F1E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of two faces, sharing exactly one edge
  void generate_non_manifold_tet_2F1E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of three tets sharing one vertex among all of them
  // and each pair sharing one edge. The three shared edges are connected in the one shared vertex. The three tets form
  // a fourth tet, which is no cell. Three of its faces belong to one of the tets each, one "face" is missing, with
  // its edges belonging to one tet each.
  void generate_non_manifold_tet_3T1V3E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of a single vertex
  void generate_non_manifold_tet_1V(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of a single edge
  void generate_non_manifold_tet_1E(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of a single face
  void generate_non_manifold_tet_1F(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of one tet and one face sharing exactly one vertex
  void generate_non_manifold_tet_1T1F1V(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of two faces, sharing exactly one vertex
  void generate_non_manifold_tet_2F1V(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of one face and one edge, sharing exactly one vertex
  void generate_non_manifold_tet_1F1E1V(TetrahedralMesh& _mesh);
  // Generate a non-manifold tetrahedral _mesh, consisting of three tets, sharing exactly one face. Only vertices and edges
  void generate_non_manifold_tet_3T1F(TetrahedralMesh& _mesh);

  // This member will be accessible in all tests
  TetrahedralMesh mesh_;
};

// simple test base for TopologicalLink and TopologicalFaceSet

class TopologicalLinkBase: public testing::Test {

protected:

    typedef OpenVolumeMesh::VertexHandle    VertexHandle;
    typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
    typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
    typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
    typedef OpenVolumeMesh::FaceHandle      FaceHandle;
    typedef OpenVolumeMesh::CellHandle      CellHandle;

    // The order of the _vertices is important, as other code relies on _vertices 0 and 4 being incident to all three tets
    void generate_triTet(TetrahedralMesh& _mesh, std::vector<VertexHandle>& _vertices);
    // generates a triangle outside another triangle, where each vertex of the inner triangle is connected to two
    // vertices of the outer triangle. This leads to the outer triangle not having a face.
    void generate_triangle_without_face(TetrahedralMesh& _mesh);
};

class TopologicalFaceSetBase: public testing::Test {

protected:

    typedef OpenVolumeMesh::TopologicalFaceSet TopologicalFaceSet;

    void generate_set_of_tet(TopologicalFaceSet& _set);

    TopologicalFaceSet set;

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

