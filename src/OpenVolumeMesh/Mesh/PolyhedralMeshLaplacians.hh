#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>

namespace OpenVolumeMesh{

template<class _polyhedral_mesh>
class PolyhedralMeshLaplacian
{
public:

    PolyhedralMeshLaplacian(_polyhedral_mesh& mesh) : mesh_(mesh) {}

    double operator[](const EdgeHandle& edge){
        return 1;
    }


    double operator[](const EdgeHandle& edge) const{
        return 1;
    }

private:

    _polyhedral_mesh& mesh_;
};

}
