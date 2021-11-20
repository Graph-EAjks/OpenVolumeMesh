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
#include <OpenVolumeMesh/Attribs/TetVolumeAttrib.hh>
#include <OpenVolumeMesh/Attribs/TriangleAreaAttrib.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>

namespace OpenVolumeMesh {

template <class GeomKernelT>
class TetHeightAttrib {
    static const inline std::string prop_name = "ovm:attrib:height";
public:
    using Scalar = typename GeomKernelT::Point::value_type;

    explicit TetHeightAttrib(GeomKernelT const& _kernel)
        : kernel_(_kernel)
        , hf_height_(_kernel.template create_private_property<double, Entity::HalfFace>(
            prop_name, std::numeric_limits<Scalar>::signaling_NaN()))
    {}

    void update(TetVolumeAttrib<GeomKernelT> const& _volume,
                TriangleAreaAttrib<GeomKernelT> const& _area)
    {
        for (auto hfh: kernel_.halffaces()) {
            const auto ch = kernel_.incident_cell(hfh);
            if (!ch.is_valid()) {
                // boundary halfface
                continue;
            }
            const auto vol = _volume[kernel_.incident_cell(hfh)];
            const auto area = _area[hfh.face_handle()];
            hf_height_[hfh] = 3 * vol / area;
        }
    }


    Scalar operator[](HalfFaceHandle _hfh) const {
        return hf_height_[_hfh];
    }

private:

    void compute_vertex_normal(const VertexHandle& _vh);

    GeomKernelT const& kernel_;

    HalfFacePropertyT<Scalar> hf_height_;

};

} // Namespace OpenVolumeMesh


