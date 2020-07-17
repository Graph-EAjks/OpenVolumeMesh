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


    /** See [Alexa 2020] Properties of Laplace Operators for Tetrahedral Meshes
        for the detail of the computations and naming
*/
    Scalar halfedge_weight(const HalfEdgeHandle& edge) const{

        Scalar weight(0);

        double alpha, beta, theta;

        auto hehf_iter = this->mesh_.hehf_iter(edge);

        /*the idea is to iterate through the halffaces around the halfedge to
         * compute the successive beta angles and then use them as alpha angles
         * for the next iteration. */
        VertexHandle xi = this->mesh_.halfedge_vertices(edge)[0];
        VertexHandle xj = this->mesh_.halfedge_vertices(edge)[1];

        VertexHandle xk = *(this->mesh_.halfface_vertices(*hehf_iter).first);

        VecT xk_xi = (this->mesh_.vertex(xi) - this->mesh_.vertex(xk)).normalized();
        VecT xk_xj = (this->mesh_.vertex(xj) - this->mesh_.vertex(xk)).normalized();

        VecT nijk = xk_xi.cross(xk_xj).normalized();

        alpha = acos(xk_xi.dot(xk_xj));

        double first_alpha = alpha;
        VecT first_nijk = nijk;

        hehf_iter++;
        int cell_count(0);

        while(hehf_iter.valid()){

            VertexHandle xl = *(this->mesh_.halfface_vertices(*hehf_iter).first);

            VecT xl_xi = (this->mesh_.vertex(xi) - this->mesh_.vertex(xl)).normalized();
            VecT xl_xj = (this->mesh_.vertex(xj) - this->mesh_.vertex(xl)).normalized();

            VecT nijl = xl_xi.cross(xl_xj).normalized();

            beta = acos(xl_xi.dot(xl_xj));
            theta = acos(nijk.dot(nijl));

            double cotan_alpha = cot(alpha);
            double cotan_beta  = cot(beta);

            weight += cot(theta) * (2. * (cotan_alpha * cotan_beta)/cos(theta) - cotan_alpha*cotan_alpha - cotan_beta*cotan_beta);

            //setting-up next iteration
            xk = xl;
            alpha = beta;
            nijk = nijl;

            hehf_iter++;
            cell_count++;
        }

        //compute the weight contribution of the last cell, which consists of the
        //last halfface and the first one
        beta = first_alpha;
        theta = acos(nijk.dot(first_nijk));

        double cotan_alpha = cot(alpha);
        double cotan_beta  = cot(beta);

        weight += cot(theta) * (2. * (cotan_alpha * cotan_beta)/cos(theta) - cotan_alpha*cotan_alpha - cotan_beta*cotan_beta);

        weight *= (this->mesh_.vertex(xi) - this->mesh_.vertex(xj)).norm() / 8.;

        return weight;
    }


private:

    double cot(Scalar x) const{
        return cos(x)/sin(x);
    }


};




}
