cmake_minimum_required (VERSION 2.6)


if("${PROJECT_NAME}" STREQUAL "")
    message("Setting project name to OpenVolumeMesh")
    project (OpenVolumeMesh)

    set(CMAKE_MACOSX_RPATH 1)

endif()

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX "d")

# include our cmake files
include (ACGCommon)

acg_get_version ()

include (ACGOutput)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# ========================================================================
# Call the subdirectories with the projects
# ========================================================================

add_subdirectory (src)

# Only add examples and doc to build target if OpenVolumeMesh
# is not built as external library only
if(${PROJECT_NAME} MATCHES "OpenVolumeMesh")
    add_subdirectory (examples)
    add_subdirectory (documentation)
else()
    set (OPENVOLUMEMESH_FOUND true PARENT_SCOPE)
    set (OPENVOLUMEMESH_LIBRARIES OpenVolumeMesh PARENT_SCOPE)
    set (OPENVOLUMEMESH_LIBRARY OpenVolumeMesh PARENT_SCOPE)
    set (OPENVOLUMEMESH_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/" PARENT_SCOPE)
    set (OPENVOLUMEMESH_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The directory where the OpenVolumeMesh libraries can be found.")
endif()

# ========================================================================
# Bundle generation (Targets exist, now configure them)
# ========================================================================

# display results
acg_print_configure_header (OPENVOLUMEMESH "OpenVolumeMesh")
