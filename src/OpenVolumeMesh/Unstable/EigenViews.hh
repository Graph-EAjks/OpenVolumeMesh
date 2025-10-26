#pragma once

#include <Eigen/Core>
#include <OpenVolumeMesh/Core/Properties/PropertyStoragePtr.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <type_traits>

namespace OpenVolumeMesh {
namespace detail {


template<typename From, typename To>
using copy_const = typename std::conditional<
    std::is_const<From>::value,
    const To,
    To
>::type;

template<typename T, typename = void>
struct eigen_view_traits {
    using value_type = T;
    static constexpr Eigen::Index dimension = 1;
    static constexpr int options = Eigen::ColMajor;

    static value_type* data(T* elem) { return elem;}
};

// Specialization for VectorT types
template<typename Scalar, int DIM>
struct eigen_view_traits<Geometry::VectorT<Scalar, DIM>> {
    using T = Geometry::VectorT<Scalar, DIM>;
    using value_type = Scalar;
    static constexpr Eigen::Index dimension = DIM;
    static constexpr int options = Eigen::RowMajor;

    static value_type* data(T* elem) { return elem->data();}
};

// Specialization for Eigen::Matrix types (column or row vectors)
template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
struct eigen_view_traits<Eigen::Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>
,typename std::enable_if<Rows_ >= 0 && Cols_ >= 0 >::type
    >
{
    using T = Eigen::Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>;
    // TODO: enable only if neither rows nor cols are -1
    using value_type = Scalar_;
    static constexpr Eigen::Index dimension = (Rows_ == 1) ? Cols_ : Rows_;
    static constexpr int options = Eigen::RowMajor;

    static value_type* data(T* elem) { return elem->data();}
};

/**
 * @brief Creates an Eigen Map view for PropertyStoragePtr<T>.
 *
 * This function provides a zero-copy view of the property storage as an Eigen vector or matrix:
 * - For scalar types (int, float, double): Returns a column vector (Dynamic x 1)
 * - For VectorT types: Returns a row-major matrix (Dynamic x DIM)
 * - For Eigen::Matrix types (vectors): Returns a row-major matrix (Dynamic x DIM)
 *
 * The constness of the returned view matches the constness of the input property.
 * Template arguments are automatically deduced from the PropertyStoragePtr type.
 *
 * @tparam PropT PropertyStoragePtr<T> or const PropertyStoragePtr<T>
 * @param prop Property storage pointer
 * @return Eigen::Map view with appropriate constness
 *
 * @note The returned view is only valid as long as the underlying property storage
 *       is not resized or destroyed.
 *
 * Example usage:
 * @code
 * // Scalar property
 * PropertyStoragePtr<double> values;
 * auto view = eigen_view(values);  // Type deduced automatically!
 * view.array() += 1.0;  // Add 1 to all values
 *
 * // VectorT property
 * PropertyStoragePtr<Vec3d> positions;
 * auto view2 = eigen_view(positions);
 * view2.col(2).array() += 1.0;  // Add 1 to all z-coordinates
 * Eigen::Vector3d mean = view2.colwise().mean();  // Compute centroid
 *
 * // Eigen vector property
 * PropertyStoragePtr<Eigen::Vector3d> normals;
 * auto view3 = eigen_view(normals);
 * view3.rowwise().normalize();  // Normalize all vectors
 * @endcode
 */
template<typename PropT>
auto eigen_view(PropT& prop) ->
        Eigen::Map<copy_const<PropT, Eigen::Matrix<
            typename eigen_view_traits<typename PropT::value_type>::value_type,
            Eigen::Dynamic,
            eigen_view_traits<typename PropT::value_type>::dimension,
            eigen_view_traits<typename PropT::value_type>::options
            >>>
{
    using ValueType = typename PropT::value_type;
    using Traits = eigen_view_traits<ValueType>;

    return {
        Traits::data(prop.data_vector().empty() ? nullptr: prop.data_vector().data()),
        Eigen::Index(prop.size()),
        Traits::dimension
    };
}
} // namespace detail
namespace Unstable {
    using OpenVolumeMesh::detail::eigen_view;
}

} // namespace OpenVolumeMesh
