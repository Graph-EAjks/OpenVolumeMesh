#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <cfloat>


namespace OpenVolumeMesh{

namespace Laplacian{



/** \brief Gives all-1 halfedge weights*/
template<class _polyhedral_mesh>
class UniformEdgeWeightEvaluator{

public:

    using Scalar = typename _polyhedral_mesh::PointT::value_type;

    static Scalar halfedge_weight(_polyhedral_mesh& mesh,
                                  const HalfEdgeHandle& edge){
        return 1.;
    }

};




/** \brief Laplacian based on the "Dual Face" of edges.
 * See [Alexa 2020] Properties of Laplace Operators for Tetrahedral Meshes. */
template<class _tetrahedral_mesh>
class DualEdgeWeightEvaluator{


public:

    using VecT   = typename _tetrahedral_mesh::PointT;
    using Scalar = typename VecT::value_type;


    /** See paper mentioned above for computations and naming detail */
    static Scalar halfedge_weight(_tetrahedral_mesh& mesh,
                                  const HalfEdgeHandle& edge){

        Scalar weight(0);

        double alpha(0), beta(0), theta(0);

        auto hehf_iter = mesh.hehf_iter(edge);

        /*the idea is to iterate through the halffaces around the halfedge to
         * compute the successive beta angles and then use them as alpha angles
         * for the next iteration. */
        VertexHandle xi = mesh.halfedge_vertices(edge)[0];
        VertexHandle xj = mesh.halfedge_vertices(edge)[1];

        VertexHandle xk;
        for(const auto& v_it: mesh.halfface_vertices(*hehf_iter)){
            if(v_it != xi && v_it != xj){
                xk = v_it;
            }
        }

        VecT xk_xi = (mesh.vertex(xi) - mesh.vertex(xk)).normalized();
        VecT xk_xj = (mesh.vertex(xj) - mesh.vertex(xk)).normalized();

        VecT nijk = xk_xi.cross(xk_xj).normalized();

        alpha = acos(xk_xi.dot(xk_xj));

        double first_alpha = alpha;
        VecT first_nijk = nijk;

        hehf_iter++;

        while(hehf_iter.valid()){

            VertexHandle xl;
            for(const auto& v_it: mesh.halfface_vertices(*hehf_iter)){
                if(v_it != xi && v_it != xj){
                    xl = v_it;
                }
            }

            VecT xl_xi = (mesh.vertex(xi) - mesh.vertex(xl)).normalized();
            VecT xl_xj = (mesh.vertex(xj) - mesh.vertex(xl)).normalized();

            VecT nijl = xl_xi.cross(xl_xj).normalized();

            beta = acos(xl_xi.dot(xl_xj));
            theta = acos(nijk.dot(nijl));

            weight += tet_weight(alpha, beta, theta);

            //setting-up next iteration
            xk = xl;
            alpha = beta;
            nijk = nijl;

            hehf_iter++;
        }

        //compute the weight contribution of the last cell, which consists of the
        //last halfface and the first one
        beta = first_alpha;
        theta = acos(nijk.dot(first_nijk));

        weight += tet_weight(alpha, beta, theta);

        //multiply by Vol(i,j) (i.e. edge length)
        weight *= (mesh.vertex(xi) - mesh.vertex(xj)).norm() / 8.;

        return weight;
    }


private:

    static Scalar cot(Scalar x){
        return cos(x)/sin(x);
    }


    static Scalar tet_weight(Scalar alpha, Scalar beta, Scalar theta) {

        if(fabs(sin(alpha)) > std::numeric_limits<Scalar>::epsilon()  &&
                fabs(sin(beta)) > std::numeric_limits<Scalar>::epsilon() &&
                fabs(sin(theta)) > std::numeric_limits<Scalar>::epsilon() &&
                fabs(cos(theta)) > std::numeric_limits<Scalar>::epsilon()){

            double cotan_alpha = cot(alpha);
            double cotan_beta  = cot(beta);

            //weight computation following the paper's formula
            return cot(theta) * (2. * (cotan_alpha * cotan_beta)/cos(theta) - cotan_alpha*cotan_alpha - cotan_beta*cotan_beta);;
        }else{
            return 0.f;
        }
    }

};


}
}

