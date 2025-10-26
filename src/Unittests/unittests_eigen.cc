#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

#ifdef __clang__
#  pragma GCC diagnostic ignored "-Weverything"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wglobal-constructors"
#  pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#  pragma GCC diagnostic ignored "-Wmissing-noreturn"
#endif

#include <gtest/gtest.h>

#include <OpenVolumeMesh/Unstable/EigenViews.hh>
#include <Eigen/Core>

using namespace OpenVolumeMesh;
using namespace Geometry;

template <typename Scalar>
class EigenTest : public testing::Test {
    public:
    static constexpr int DIM = 3;
    using Vec11T = VectorT<Scalar, DIM>;
    using EigenVec = Eigen::Vector<Scalar, DIM>;
    using Mesh = TopologyKernel;
    Mesh mesh;
    VH v0 = mesh.add_vertex();
    VH v1 = mesh.add_vertex();
};

// Define the types to test
using VectorTypes = testing::Types<float, double, int>;
TYPED_TEST_SUITE(EigenTest, VectorTypes);

TYPED_TEST(EigenTest, PositionView) {
    using Scalar = TypeParam;
    using Vec11T = VectorT<Scalar, 3>;
    using EigenVec = Eigen::Vector<Scalar, 3>;
    using Mesh = GeometryKernel<Vec11T, TopologyKernel>;
    Mesh mesh;
    auto v0 = mesh.add_vertex({1,2,3});
    auto v1 = mesh.add_vertex({4,5,6});

    auto X = Unstable::eigen_view(mesh.vertex_positions());

    ASSERT_EQ(X.rows(), mesh.n_vertices());
    ASSERT_EQ(X.cols(), 3);
    EXPECT_EQ(X.row(0), (EigenVec() << 1,2,3).finished().transpose());
    EXPECT_EQ(X.row(1), (EigenVec() << 4,5,6).finished().transpose());

    X *=2;
    EXPECT_EQ(mesh.vertex(v0), Vec11T(2,4,6));
}

TYPED_TEST(EigenTest, ViewScalarProps) {
    using Scalar = TypeParam;
    using Mesh = TopologyKernel;
    auto prop = this->mesh.template create_private_vertex_property<Scalar>("flt");
    prop[VH(0)] = 23;
    prop[VH(1)] = 42;

    auto F = Unstable::eigen_view(prop);
    EXPECT_EQ(prop[VH(0)], F[0]);
    EXPECT_EQ(prop[VH(1)], F[1]);
    F *= 2;
    EXPECT_EQ(prop[VH(0)], 46);
    EXPECT_EQ(prop[VH(0)], F[0]);
}

TYPED_TEST(EigenTest, EigenDynamicVectorProps) {
    using Vec = Eigen::Vector<TypeParam, Eigen::Dynamic>;
    auto prop_eigen_dyn_vec = this->mesh.template create_private_vertex_property<Vec>();
    prop_eigen_dyn_vec[VH(0)].setZero(4);
    prop_eigen_dyn_vec[VH(1)].setOnes(2);
    auto dv = Unstable::eigen_view(prop_eigen_dyn_vec);
    ASSERT_EQ(dv.rows(), this->mesh.n_vertices());
    ASSERT_EQ(dv.cols(), 1); // an eigen vector of vectorXd..
    ASSERT_EQ(dv[0].cols(), 1);
    ASSERT_EQ(dv[0].rows(), 4);
    ASSERT_EQ(dv[1].cols(), 1);
    ASSERT_EQ(dv[1].rows(), 2);
}
TYPED_TEST(EigenTest, EigenVectorRowVectorProps) {
    using Scalar = TypeParam;
    auto prop_eigen_row_vec = this->mesh.template create_private_vertex_property<Eigen::RowVector<Scalar, 4>>();
    prop_eigen_row_vec[VH(0)] = {23, 42, 0, 0};
    prop_eigen_row_vec[VH(1)] = {10, 50, 0, 0};
    auto rv = Unstable::eigen_view(prop_eigen_row_vec);

    ASSERT_EQ(rv.rows(), this->mesh.n_vertices());
    ASSERT_EQ(rv.cols(), 4);
    ASSERT_EQ(rv(0, 0), 23);
    ASSERT_EQ(rv(0, 1), 42);
    ASSERT_EQ(rv(1, 0), 10);
    ASSERT_EQ(rv(1, 1), 50);

}
TYPED_TEST(EigenTest, EigenVectorColVectorProps) {
    using Scalar = TypeParam;
    auto prop_eigen_col_vec = this->mesh.template create_private_vertex_property<Eigen::Vector<Scalar, 4>>();
    prop_eigen_col_vec[VH(0)] = {23, 42, 0, 0};
    prop_eigen_col_vec[VH(1)] = {10, 50, 0, 0};
    auto cv = Unstable::eigen_view(prop_eigen_col_vec);

    ASSERT_EQ(cv.rows(), this->mesh.n_vertices());
    ASSERT_EQ(cv.cols(), 4);
    ASSERT_EQ(cv(0, 0), 23);
    ASSERT_EQ(cv(0, 1), 42);
    ASSERT_EQ(cv(1, 0), 10);
    ASSERT_EQ(cv(1, 1), 50);



}

