#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>

namespace OpenVolumeMesh{


/**\brief base class for all Laplacians
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
    BaseLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh){}


    virtual double halfedge_weight(const HalfEdgeHandle& edge) const = 0;

protected:

    _polyhedral_mesh& mesh_;

};


template<class _polyhedral_mesh>
class UniformLaplacian : public BaseLaplacian<_polyhedral_mesh>
{
public:

    UniformLaplacian(_polyhedral_mesh& mesh) : BaseLaplacian<_polyhedral_mesh>(mesh) {}

    double halfedge_weight(const HalfEdgeHandle& edge) const{
        return 1;
    }

};


template< template<class> class _base_laplacian, class _polyhedral_mesh>
class Laplacian : public _base_laplacian<_polyhedral_mesh>{

public:

    using BaseLaplacian = _base_laplacian<_polyhedral_mesh>;

    Laplacian(_polyhedral_mesh& mesh) : BaseLaplacian(mesh) {}


    double operator[](const HalfEdgeHandle& half_edge) const {
        return BaseLaplacian::halfedge_weight(half_edge);
    }

    Vec3d operator[](const VertexHandle& vertex) const{

        Vec3d weighted_average({0,0,0});
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



}
