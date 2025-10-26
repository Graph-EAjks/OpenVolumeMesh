#pragma once
/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/


#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Core/AddSmartIterators.hh>

namespace OpenVolumeMesh {


class OVM_EXPORT PolyhedralMeshTopologyKernel
    : public AddSmartIterators<PolyhedralMeshTopologyKernel, TopologyKernel>
{
};


/*
 * Predefines for most common mesh types
 */
typedef GeometryKernel<Geometry::Vec2i, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2i; typedef GeometryKernel<Geometry::Vec2ui, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2ui;
typedef GeometryKernel<Geometry::Vec2f, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2f;
typedef GeometryKernel<Geometry::Vec2d, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2d;
typedef GeometryKernel<Geometry::Vec2c, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2c;
typedef GeometryKernel<Geometry::Vec2uc, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV2uc;
typedef GeometryKernel<Geometry::Vec3i, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3i;
typedef GeometryKernel<Geometry::Vec3ui, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3ui;
typedef GeometryKernel<Geometry::Vec3f, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3f;
typedef GeometryKernel<Geometry::Vec3d, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3d;
typedef GeometryKernel<Geometry::Vec3c, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3c;
typedef GeometryKernel<Geometry::Vec3uc, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV3uc;
typedef GeometryKernel<Geometry::Vec4i, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4i;
typedef GeometryKernel<Geometry::Vec4ui, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4ui;
typedef GeometryKernel<Geometry::Vec4f, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4f;
typedef GeometryKernel<Geometry::Vec4d, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4d;
typedef GeometryKernel<Geometry::Vec4c, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4c;
typedef GeometryKernel<Geometry::Vec4uc, PolyhedralMeshTopologyKernel> GeometricPolyhedralMeshV4uc;

using TopologicPolyhedralMesh = PolyhedralMeshTopologyKernel;

} // Namespace OpenVolumeMesh

