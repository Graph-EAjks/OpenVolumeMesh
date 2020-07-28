#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <cfloat>


#pragma message("replace with static const")
#define ZERO_SIN_THETA_THRESHOLD DBL_EPSILON

namespace OpenVolumeMesh{


/** \brief generic template for all Laplacians
 *  _base_laplacian should extend the BaseLaplacian template (see below)
*/
template< template<class> class _base_vertex_laplacian, class _polyhedral_mesh>
class VertexLaplacian : public _base_vertex_laplacian<_polyhedral_mesh>{

public:

    using BaseVertexLaplacian = _base_vertex_laplacian<_polyhedral_mesh>;
    using VecT   = typename _polyhedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    VertexLaplacian(_polyhedral_mesh& mesh) : BaseVertexLaplacian(mesh) {}


    Scalar operator[](const HalfEdgeHandle& half_edge) const {
        return BaseVertexLaplacian::halfedge_weight(half_edge);
    }


    VecT operator[](const VertexHandle& vertex) const{

        VecT weighted_sum = {0,0,0};
        double weight_sum(0);

        auto voh_it = this->mesh_.voh_iter(vertex);
        while(voh_it.valid()){

            double he_weight = BaseVertexLaplacian::halfedge_weight(*voh_it);

            weighted_sum += he_weight * this->mesh_.vertex(this->mesh_.to_vertex_handle(*voh_it));
            weight_sum += he_weight;

            if(weight_sum != weight_sum){
                exit(1);
            }

            voh_it++;
        }

        return weighted_sum / weight_sum;
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
class BaseVertexLaplacian{

public:

    using Scalar = typename _polyhedral_mesh::PointT::value_type;

    BaseVertexLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh){}


    virtual Scalar halfedge_weight(const HalfEdgeHandle& edge) const = 0;

protected:

    _polyhedral_mesh& mesh_;

};



/** \brief Gives all-1 halfedge weights*/
template<class _polyhedral_mesh>
class UniformVertexLaplacian : public BaseVertexLaplacian<_polyhedral_mesh>{

public:


    UniformVertexLaplacian(_polyhedral_mesh& mesh) : BaseVertexLaplacian<_polyhedral_mesh>(mesh) {}

    double halfedge_weight(const HalfEdgeHandle& edge) const{
        return 1.;
    }

};





/** \brief Laplacian based on the "Dual Face" of edges.
 * See [Alexa 2020] Properties of Laplace Operators for Tetrahedral Meshes. */
template<class _tetrahedral_mesh>
class DualLaplacian : public BaseVertexLaplacian<_tetrahedral_mesh>{


public:

    using VecT   = typename _tetrahedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    DualLaplacian(_tetrahedral_mesh& mesh) : BaseVertexLaplacian<_tetrahedral_mesh>(mesh) {}


    /** See paper mentioned above for computations and naming detail */
    Scalar halfedge_weight(const HalfEdgeHandle& edge) const{

        Scalar weight(0);

        double alpha, beta, theta;

        auto hehf_iter = this->mesh_.hehf_iter(edge);

        /*the idea is to iterate through the halffaces around the halfedge to
         * compute the successive beta angles and then use them as alpha angles
         * for the next iteration. */
        VertexHandle xi = this->mesh_.halfedge_vertices(edge)[0];
        VertexHandle xj = this->mesh_.halfedge_vertices(edge)[1];

        VertexHandle xk;
        for(auto v_it = this->mesh_.halfface_vertices(*hehf_iter).first;
            v_it != this->mesh_.halfface_vertices(*hehf_iter).second;
            v_it++){
            if(*v_it != xi && *v_it != xj){
                xk = *v_it;
            }
        }

        VecT xk_xi = (this->mesh_.vertex(xi) - this->mesh_.vertex(xk)).normalized();
        VecT xk_xj = (this->mesh_.vertex(xj) - this->mesh_.vertex(xk)).normalized();

        VecT nijk = xk_xi.cross(xk_xj).normalized();

        alpha = acos(xk_xi.dot(xk_xj));

        double first_alpha = alpha;
        VecT first_nijk = nijk;

        hehf_iter++;
        int cell_count(0);

        while(hehf_iter.valid()){

            VertexHandle xl;
            for(auto v_it = this->mesh_.halfface_vertices(*hehf_iter).first;
                v_it != this->mesh_.halfface_vertices(*hehf_iter).second;
                v_it++){
                if(*v_it != xi && *v_it != xj){
                    xl = *v_it;
                }
            }

            VecT xl_xi = (this->mesh_.vertex(xi) - this->mesh_.vertex(xl)).normalized();
            VecT xl_xj = (this->mesh_.vertex(xj) - this->mesh_.vertex(xl)).normalized();

            VecT nijl = xl_xi.cross(xl_xj).normalized();

            beta = acos(xl_xi.dot(xl_xj));
            theta = acos(nijk.dot(nijl));

            weight += tet_weight(alpha, beta, theta);

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

        weight += tet_weight(alpha, beta, theta);

        //multiply by Vol(i,j) (edge length)
        weight *= (this->mesh_.vertex(xi) - this->mesh_.vertex(xj)).norm() / 8.;

        return weight;
    }


private:

    double cot(Scalar x) const{
        return cos(x)/sin(x);
    }


    double tet_weight(double alpha, double beta, double theta) const {
        if(abs(sin(alpha)) > ZERO_SIN_THETA_THRESHOLD  &&
                abs(sin(beta)) > ZERO_SIN_THETA_THRESHOLD &&
                abs(sin(theta)) > ZERO_SIN_THETA_THRESHOLD &&
                abs(cos(theta)) > ZERO_SIN_THETA_THRESHOLD){

            double cotan_alpha = cot(alpha);
            double cotan_beta  = cot(beta);

            //weight computation following the paper's formula
            return cot(theta) * (2. * (cotan_alpha * cotan_beta)/cos(theta) - cotan_alpha*cotan_alpha - cotan_beta*cotan_beta);;
        }else{
            return 0.f;
        }
    }

};




template< template<class> class _base_laplacian, class _polyhedral_mesh>
class PrecomputedLaplacian : public _base_laplacian<_polyhedral_mesh>{

public:

    using BaseLaplacian = _base_laplacian<_polyhedral_mesh>;
    using VecT   = typename _polyhedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    PrecomputedLaplacian(_polyhedral_mesh& mesh) :
        BaseLaplacian(mesh),
        vertex_laplacians_(mesh. template request_vertex_property<VecT>("laplacians")),
        edge_weights_(mesh. template request_halfedge_property<Scalar>("laplacian weights")){
        //pre-computations


    }

    Scalar operator[](const HalfEdgeHandle& edge) const{
        return 0.f;
    }

    VecT operator[](const VertexHandle& vertex) const{
        return {0,0,0};

    }


private:

    VertexPropertyT<VecT> vertex_laplacians_;
    HalfEdgePropertyT<Scalar> edge_weights_;
};



}
