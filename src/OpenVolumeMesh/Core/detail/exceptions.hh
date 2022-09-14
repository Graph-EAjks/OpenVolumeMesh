#pragma once
#include <OpenVolumeMesh/Config/Export.hh>
#include <stdexcept>

namespace OpenVolumeMesh::Core::detail {

    class OVM_EXPORT core_error : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

class OVM_EXPOR invalid_topology_error : public core_error
{
using io_error::io_error;
};


} // namespace OpenVolumeMesh::Core::detail
