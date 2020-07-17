#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

namespace OpenVolumeMesh{


/** \brief generic template for all Laplacians
 *  _base_laplacian should extend the BaseLaplacian template (see below)
*/
template< template<class> class _base_laplacian, class _polyhedral_mesh>
class Laplacian : public _base_laplacian<_polyhedral_mesh>{

public:

    using BaseLaplacian = _base_laplacian<_polyhedral_mesh>;
    using VecT   = typename _polyhedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    Laplacian(_polyhedral_mesh& mesh) : BaseLaplacian(mesh) {}


    Scalar operator[](const HalfEdgeHandle& half_edge) const {
        return BaseLaplacian::halfedge_weight(half_edge);
    }

    VecT operator[](const VertexHandle& vertex) const{

        VecT weighted_average({0,0,0});
        double weight_sum(0);

        auto voh_it = this->mesh_.voh_iter(vertex);
        while(voh_it.valid()){

            double he_weight = BaseLaplacian::halfedge_weight(*voh_it);

            weighted_average += he_weight * this->mesh_.vertex(this->mesh_.to_vertex_handle(*voh_it));
            weight_sum += he_weight;

            voh_it++;
        }
        return weighted_average / weight_sum;
    }


private:


};

/** \brief base class for all Laplacians
 *
 * To create your own Laplacian, make it inherit from this base class and
 * implement the 'half_edge_weight' function with the desired weight computation
 *
 * You can then use it as template argument for the 'Laplacian' template.
 * See UniformLaplacian for an example
*/
template<class _polyhedral_mesh>
class BaseLaplacian{

public:

    using Scalar = typename _polyhedral_mesh::PointT::value_type;

    BaseLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh){}


    virtual Scalar halfedge_weight(const HalfEdgeHandle& edge) const = 0;

protected:

    _polyhedral_mesh& mesh_;

};


template<class _polyhedral_mesh>
class UniformLaplacian : public BaseLaplacian<_polyhedral_mesh>{

public:


    UniformLaplacian(_polyhedral_mesh& mesh) : BaseLaplacian<_polyhedral_mesh>(mesh) {}

    double halfedge_weight(const HalfEdgeHandle& edge) const{
        return 1;
    }

};




template<class _tetrahedral_mesh>
class DualLaplacian : public BaseLaplacian<_tetrahedral_mesh>{

    //static_assert(_tetrahedral_mesh::)

public:

    using VecT   = typename _tetrahedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    DualLaplacian(_tetrahedral_mesh& mesh) : BaseLaplacian<_tetrahedral_mesh>(mesh) {}

    Scalar halfedge_weight(const HalfEdgeHandle& edge) const{


        return 1;
    }


private:

    VecT circumcenter(const VecT& a,
                      const VecT& b,
                      const VecT& c){

        VecT cc;

        const double l[3]{
            (b - c).squaredNorm(),
            (a - c).squaredNorm(),
            (a - b).squaredNorm()
        };

        const double ba[3]{l[0] * (l[1] + l[2] - l[0]), l[1] * (l[2] + l[0] - l[1]), l[2] * (l[0] + l[1] - l[2])};
        const double sum = ba[0] + ba[1] + ba[2];

        cc = (ba[0] / sum) * a + (ba[1] / sum) * b + (ba[2] / sum) * c;
    }

};




}
