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

#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <type_traits>

namespace OpenVolumeMesh {

ResourceManager::ResourceManager(const ResourceManager &other)
{
    copyAllPropertiesFrom(other);
    *this = other;
}

ResourceManager& ResourceManager::operator=(ResourceManager &&other)
{
    if (this == &other) return *this;
    moveAllPropertiesFrom(std::move(other));
    return *this;
}


ResourceManager& ResourceManager::operator=(const ResourceManager &other)
{
    if (this == &other) return *this;
    copyAllPropertiesFrom(other);
    return *this;
}

void ResourceManager::copyAllPropertiesFrom(ResourceManager const&src)
{
    PerEntityStorageTrackers result;

    for_each_entity([&](auto entity_tag) {
        using ET = decltype(entity_tag);
        clone_props(
                    src.storage_tracker<ET>(),
                    result.get<ET>(),
                    persistent_props_.get<ET>());
        merge_props(
                    storage_tracker<ET>(),
                    result.get<ET>(),
                    persistent_props_.get<ET>(),
                    src.n<ET>());
    });
    storage_trackers_ = std::move(result);
}
void ResourceManager::moveAllPropertiesFrom(ResourceManager &&src)
{

    PerEntityStorageTrackers result = std::move(src.storage_trackers_);

    for_each_entity([&](auto entity_tag) {
        using ET = decltype(entity_tag);
        persistent_props_.get<ET>().merge(src.persistent_props_.get<ET>());
        assert(src.persistent_props_.get<ET>().empty()); // there should be no common props!
        merge_props(
                    storage_tracker<ET>(),
                    result.get<ET>(),
                    persistent_props_.get<ET>(),
                    src.n<ET>());
    });
    storage_trackers_ = std::move(result);
}


detail::Tracker<PropertyStorageBase> &
ResourceManager::storage_tracker(EntityType type) const
{
    return storage_trackers_.get(type);
}

void ResourceManager::resize_vprops(size_t _nv) {
    resize_props<Entity::Vertex>(_nv);
}

void ResourceManager::resize_eprops(size_t _ne) {
    resize_props<Entity::Edge>(_ne);
    resize_props<Entity::HalfEdge>(2 * _ne);
}

void ResourceManager::resize_fprops(size_t _nf) {
    resize_props<Entity::Face>(_nf);
    resize_props<Entity::HalfFace>(2 * _nf);
}

void ResourceManager::resize_cprops(size_t _nc) {
    resize_props<Entity::Cell>(_nc);
}

void ResourceManager::reserve_vprops(size_t _n) {
    reserve_props<Entity::Vertex>(_n);
}
void ResourceManager::reserve_eprops(size_t _n) {
    reserve_props<Entity::Edge>(_n);
    reserve_props<Entity::HalfEdge>(2 * _n);
}
void ResourceManager::reserve_fprops(size_t _n) {
    reserve_props<Entity::Face>(_n);
    reserve_props<Entity::HalfFace>(2 * _n);
}
void ResourceManager::reserve_cprops(size_t _n) {
    reserve_props<Entity::Cell>(_n);
}


void ResourceManager::vertex_deleted(const VertexHandle& _h) {
    entity_deleted(_h);
}

void ResourceManager::edge_deleted(const EdgeHandle& _h) {
    entity_deleted(_h);
    entity_deleted(_h.halfedge_handle(1));
    entity_deleted(_h.halfedge_handle(0));
}

void ResourceManager::face_deleted(const FaceHandle& _h)
{
    entity_deleted(_h);
    entity_deleted(_h.halfface_handle(1));
    entity_deleted(_h.halfface_handle(0));
}

void ResourceManager::cell_deleted(const CellHandle& _h) {
    entity_deleted(_h);
}

void ResourceManager::clear_all_props()
{
    for_each_entity([this](auto entity_tag){ clear_props<decltype(entity_tag)>();});
}


template<> size_t OVM_EXPORT ResourceManager::n<Entity::Vertex>()   const { return n_vertices(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Edge>()     const { return n_edges(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::HalfEdge>() const { return n_halfedges(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Face>()     const { return n_faces(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::HalfFace>() const { return n_halffaces(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Cell>()     const { return n_cells(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Mesh>()     const { return 1; }

/// clone shared properties from src tracker into dst tracker,
/// add to our persistent set if they were persistent.
void ResourceManager::clone_props(
        StorageTracker const&src,
        StorageTracker &dst,
        PersistentProperties &persistent)
{
    for (const auto &srcprop: src) {
        if (!srcprop->shared()) {
            // it does not make sense to clone private props,
            // noone could access them
            continue;
        }
        auto dstprop = srcprop->clone();
        dstprop->set_tracker(&dst);
        if (srcprop->persistent()) {
            persistent.insert(dstprop->shared_from_this());
        }
    }
}

/// Merge existing props into the result for the moved/copied-into object (*this).
///
/// if a src property is
///     - not in dst: resize, move to dst
///     - in dst: - update PropertyPtr instances of src to point to the replacement prop
///               - remove old prop from persistent
void ResourceManager::merge_props(
        StorageTracker &src,
        StorageTracker &dst,
        PersistentProperties &persistent,
        size_t n_elem)
{
    std::vector<PropertyStorageBase*> to_be_moved;

    for (auto* srcprop: src) {
        if (!srcprop->shared()) {
            // private properties cannot match with props already in dst
            srcprop->resize(n_elem);
            to_be_moved.push_back(srcprop);
        }

        bool found = false;
        for (auto* dstprop: dst)
        {
            if (!dstprop->shared()
                    && dstprop->name() == srcprop->name()
                    && dstprop->internal_type_name() == srcprop->internal_type_name())
            {
                // found a correspondence!
                // TODO: update PropertyPtr instances here!

                persistent.erase(srcprop->shared_from_this());
                found = true;
                break;
            }
        }
        if (!found) {
            srcprop->resize(n_elem);
            to_be_moved.push_back(srcprop);
        }
    }
    for (auto& prop: to_be_moved) {
        prop->set_tracker(&dst);
    }
}



} // Namespace OpenVolumeMesh
