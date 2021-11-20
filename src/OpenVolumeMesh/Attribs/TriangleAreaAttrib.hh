#pragma once
/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2021 by Computer Graphics Group, RWTH Aachen         *
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


#include <cassert>

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Geometry/Area.hh>

namespace OpenVolumeMesh {

template <class GeomKernelT>
class TriangleAreaAttrib {
    static const inline std::string prop_name = "ovm:attrib:triangle_area";
public:
    using Scalar = typename GeomKernelT::Point::value_type;

    explicit TriangleAreaAttrib(GeomKernelT const& _kernel)
        : kernel_(_kernel)
        , face_area_(_kernel.template create_private_property<double, Entity::Face>(
            prop_name, std::numeric_limits<Scalar>::signaling_NaN()))
    {}


    void update()
    {
        for (auto fh: kernel_.faces()) {
            face_area_[fh] = compute_triangle_area(kernel_, fh);
        }
    }


    Scalar operator[](FaceHandle _fh) const {
        return face_area_[_fh];
    }

private:
    GeomKernelT const& kernel_;
    FacePropertyT<Scalar> face_area_;

};

} // Namespace OpenVolumeMesh



