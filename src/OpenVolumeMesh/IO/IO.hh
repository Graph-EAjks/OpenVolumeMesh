#pragma once

#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <filesystem>
#include <stdexcept>

namespace OpenVolumeMesh::IO {

/// Read a mesh in ovm or ovmb file format depending on the file ending.
template<typename MeshT>
[[deprecated("Prefer read_mesh()")]]
bool read_file(std::string const&_filename, MeshT &_mesh,
               bool _topo_check = true, bool _bottom_up_incidences = true)
{
    ReadOptions options;
    options.topology_check = _topo_check;
    options.bottom_up_incidences = _bottom_up_incidences;
    auto res = read_mesh(_filename, _mesh, options);
    return res == ReadResult::Ok;
}

template<typename MeshT>
ReadResult read_mesh(
        std::string const&_filename,
        MeshT &_mesh,
        ReadOptions const&options)
{
    const auto path = std::filesystem::path(_filename);
    const std::string ext = path.extension();

    if (ext == "ovmb") {
        return ovmb_read(_filename.c_str(), _mesh, options);
    } else if (ext == "ovm") {
        FileManager file_manager;
        auto ok = file_manager.readFile(_filename, _mesh,
                options.topology_check,
                options.bottom_up_incidences);
        return ok ? ReadResult::Ok : ReadResult::OtherError;
    } else {
        return ReadResult::UnknownExtension;
    }
}

template<typename MeshT>
WriteResult write_mesh(
        std::string const&_filename,
        MeshT const &_mesh,
        ReadOptions const&options)
{
    const auto path = std::filesystem::path(_filename);
    const std::string ext = path.extension();

    if (ext == "ovmb") {
        return ovmb_write(_filename.c_str(), _mesh, options);
    } else if (ext == "ovm") {
        FileManager file_manager;
        auto ok = file_manager.writeFile(_filename, _mesh);
        return ok ? WriteResult::Ok : WriteResult::Error;
    } else {
        return WriteResult::UnknownExtension;
    }
}

} // namespace OpenVolumeMesh::IO
