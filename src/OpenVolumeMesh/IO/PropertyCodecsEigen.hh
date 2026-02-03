#pragma once

#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/PropertyCodecsT_impl.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>

#include <Eigen/Core>


namespace OpenVolumeMesh::Geometry {
template<typename _Scalar, int _Cols>
struct vector_dim<Eigen::Matrix<_Scalar, 1, _Cols>> {
    constexpr static int dim = _Cols;
};
template<typename _Scalar, int _Rows>
struct vector_dim<Eigen::Matrix<_Scalar, _Rows, 1>> {
    constexpr static int dim = _Rows;
};
}

namespace OpenVolumeMesh::IO::Codecs {

template<typename _Scalar, int _Rows, int _Cols>
struct EigenDenseFixedMatrix
{
    // TODO: replace with register_matrixlike
    static_assert(_Rows > 0);
    static_assert(_Cols > 0);
    using T = Eigen::Matrix<_Scalar, _Rows, _Cols>;

    static void encode(detail::Encoder &enc, const T &val) {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                enc.write(val(r, c));
            }
        }
    }

    static void decode(detail::Decoder &reader, T &val) {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                reader.read(val(r, c));
            }
        }
    }
};

} // namespace OpenVolumeMesh::IO::Codecs

namespace OpenVolumeMesh::IO {

template<typename _Scalar, int _Rows, int _Cols>
inline void register_eigen_fixed_matrix_codec(PropertyCodecs &_codecs, std::string _name)
{
    using namespace Codecs;
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<_Scalar, _Rows, _Cols>>>(
            std::move(_name),
            Eigen::Matrix<_Scalar, _Rows, _Cols>::Zero());
}
inline void register_eigen_codecs(PropertyCodecs &_codecs)
{
    using namespace Codecs;
    register_eigen_fixed_matrix_codec<double, 2, 1>(_codecs, "2d");
    register_eigen_fixed_matrix_codec<double, 3, 1>(_codecs, "3d");
    register_eigen_fixed_matrix_codec<double, 4, 1>(_codecs, "4d");
    register_eigen_fixed_matrix_codec<double, 9, 1>(_codecs, "9d");

    register_eigen_fixed_matrix_codec<float, 2, 1>(_codecs, "2f");
    register_eigen_fixed_matrix_codec<float, 3, 1>(_codecs, "3f");
    register_eigen_fixed_matrix_codec<float, 4, 1>(_codecs, "4f");
    register_eigen_fixed_matrix_codec<float, 9, 1>(_codecs, "9f");

    // matrices stored as column-major:

    register_eigen_fixed_matrix_codec<double, 2, 2>(_codecs, "2x2d");
    register_eigen_fixed_matrix_codec<double, 3, 3>(_codecs, "3x3d");
    register_eigen_fixed_matrix_codec<double, 4, 4>(_codecs, "4x4d");

    register_eigen_fixed_matrix_codec<float, 2, 2>(_codecs, "2x2f");
    register_eigen_fixed_matrix_codec<float, 3, 3>(_codecs, "3x3f");
    register_eigen_fixed_matrix_codec<float, 4, 4>(_codecs, "4x4f");

    // TODO: dynamic (dense/sparse) types
}

} // namespace OpenVolumeMesh::IO
