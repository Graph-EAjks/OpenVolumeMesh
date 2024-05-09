#pragma once

namespace OpenVolumeMesh::IO {

struct WriteOptions
{
    enum class TopologyType {
        AutoDetect,
        Polyhedral,
        Tetrahedral,
        Hexahedral,
    } topology_type = TopologyType::AutoDetect;
};

} // namespace OpenVolumeMesh::IO
