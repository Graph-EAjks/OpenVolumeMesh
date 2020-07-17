#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>

namespace OpenVolumeMesh{

template<class _polyhedral_mesh>
class PolyhedralMeshLaplacian
{
public:

    PolyhedralMeshLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh) {}

    double operator[](const EdgeHandle& edge) const{
        return 1;
    }

    double operator[](const HalfEdgeHandle& half_edge) const {
        return 1;
    }



    Vec3d operator[](const VertexHandle& vertex) const{

        Vec3d weighted_average({0,0,0});
        double weight_sum(0);

        auto voh_it = mesh_.voh_iter(vertex);
        while(voh_it.valid()){

            double he_weight = (*this)[*voh_it];

            weighted_average += he_weight * mesh_.vertex(mesh_.to_vertex_handle(*voh_it));
            weight_sum += he_weight;

            voh_it++;
        }
        return weighted_average / weight_sum;
    }


private:

    _polyhedral_mesh& mesh_;
};

}
