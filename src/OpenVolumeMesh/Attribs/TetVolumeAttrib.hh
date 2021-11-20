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
#include <OpenVolumeMesh/Geometry/Volume.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>

namespace OpenVolumeMesh {

template <class GeomKernelT>
class TetVolumeAttrib {
public:
    using Scalar = typename GeomKernelT::Point::value_type;

    explicit TetVolumeAttrib(GeomKernelT& _kernel)
        : kernel_(_kernel)
        , cell_volume_(&_kernel, "height", std::numeric_limits<Scalar>::signaling_NaN())
    {}

    void update()
    {
        for (auto ch: kernel_.cells()) {
            cell_volume_[ch] = compute_tet_volume(kernel_, ch);
        }
    }

    Scalar operator[](CellHandle _ch) const {
        return cell_volume_[_ch];
    }

private:
    GeomKernelT const& kernel_;
    CellPropertyT<Scalar> cell_volume_;
};

} // Namespace OpenVolumeMesh


