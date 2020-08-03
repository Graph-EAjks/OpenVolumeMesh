#pragma once

#include <OpenVolumeMesh/DifferentialOperators/EdgeWeightEvaluators.hh>

namespace OpenVolumeMesh{




/** \brief generic template for all Laplacians
 *  _base_laplacian should extend the BaseLaplacian template (see below) */
template<template <class> class _edge_weight_evaluator, class _polyhedral_mesh>
class VertexLaplacian {

public:

    using VecT   = typename _polyhedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    VertexLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh) {}


    Scalar operator[](const HalfEdgeHandle& half_edge) const {
        return _edge_weight_evaluator<_polyhedral_mesh>::halfedge_weight(mesh_, half_edge);
    }


    Scalar operator[](const VertexHandle& vertex) const{

        double weight_sum(0);

        auto voh_it = this->mesh_.voh_iter(vertex);
        while(voh_it.valid()){

            weight_sum += _edge_weight_evaluator<_polyhedral_mesh>::halfedge_weight(mesh_, *voh_it);

            voh_it++;
        }

        return weight_sum;
    }


protected:

    _polyhedral_mesh& mesh_;

};




/** \brief this template is equivalent to its _base_laplacian template argument.
 * The only difference being that the weights and laplacians are computed at construction
 * and can thus be accessed fast */
template< template<class> class _edge_weight_evaluator, class _polyhedral_mesh>
class PrecomputedVertexLaplacian : public VertexLaplacian<_edge_weight_evaluator, _polyhedral_mesh>{

public:

    using EdgeWeightEvaluator = _edge_weight_evaluator<_polyhedral_mesh>;
    using BaseLaplacian = VertexLaplacian<_edge_weight_evaluator, _polyhedral_mesh>;

    using VecT   = typename _polyhedral_mesh::PointT;
    using Scalar = typename VecT::value_type;

    PrecomputedVertexLaplacian(_polyhedral_mesh& mesh) :
        BaseLaplacian(mesh),
        vertex_weights_(mesh. template request_vertex_property<Scalar>("laplacians")),
        edge_weights_(mesh. template request_edge_property<Scalar>("laplacian weights")){

        //pre-computations

        //per edge
        for(const auto& edge: this->mesh_.edges()){
            edge_weights_[edge] = EdgeWeightEvaluator::halfedge_weight(this->mesh_,
                                                                       this->mesh_.halfedge_handle(edge, 0));
        }

        //and per vertex
        for(const auto& vertex: this->mesh_.vertices()){
            Scalar weight_sum = 0;

            auto voh = this->mesh_.voh_iter(vertex);
            while(voh.valid()){
                weight_sum += edge_weights_[this->mesh_.edge_handle(*voh)];;
                voh++;
            }

            vertex_weights_[vertex] = weight_sum;
        }

    }

    Scalar operator[](const HalfEdgeHandle& edge) const{
        return edge_weights_[this->mesh_.edge_handle(edge)];
    }

    Scalar operator[](const VertexHandle& vertex) const{
        return vertex_weights_[vertex];
    }


private:

    VertexPropertyT<Scalar> vertex_weights_;
    EdgePropertyT<Scalar>   edge_weights_;
};

/* Laplacian aliases */
template<class _polyhedral_mesh>
using UniformVertexLaplacian = VertexLaplacian<Laplacian::UniformEdgeWeightEvaluator, _polyhedral_mesh>;
template<class _polyhedral_mesh>
using DualVertexLaplacian = VertexLaplacian<Laplacian::DualEdgeWeightEvaluator, _polyhedral_mesh>;

/* Pre-computed Laplacian aliases */
template<class _polyhedral_mesh>
using PrecomputedUniformVertexLaplacian = PrecomputedVertexLaplacian<Laplacian::UniformEdgeWeightEvaluator, _polyhedral_mesh>;
template<class _polyhedral_mesh>
using PrecomputedDualVertexLaplacian = PrecomputedVertexLaplacian<Laplacian::DualEdgeWeightEvaluator, _polyhedral_mesh>;

}
