cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
cmake_policy(VERSION 3.24.0)

if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    find_package(OpenVolumeMesh REQUIRED)
endif()

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

add_executable(ovm2usd ovm2usd/ovm2usd.cc)

if (NOT TARGET tinyusdz::tinyusdz_static)
    FetchContent_Declare(tinyusdz
        GIT_REPOSITORY https://github.com/lighttransport/tinyusdz
        GIT_TAG ba6b79d163ebcb837cb423628ea612cdf02a82a7 # dev branch 2025-01-26
    )
    FetchContent_MakeAvailable(tinyusdz)

    target_include_directories(ovm2usd PRIVATE "$<BUILD_INTERFACE:${tinyusdz_SOURCE_DIR}/src>")
endif()


target_link_libraries(ovm2usd
    PRIVATE
    OpenVolumeMesh::OpenVolumeMesh
    tinyusdz::tinyusdz_static
    )


#[[=
if (NOT TARGET happly::happly)
    FetchContent_Declare(happly
        GIT_REPOSITORY https://github.com/nmwsharp/happly
        GIT_TAG 8a606309daaa680eee495c8279feb0b704148f4a # master branch, Feb 7, 2024
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
    )
    FetchContent_MakeAvailable(happly)
    add_library(happly::happly INTERFACE IMPORTED)
    target_include_directories(happly::happly INTERFACE "$<BUILD_INTERFACE:${happly_SOURCE_DIR}>")
    message("Downloaded happy to ${happly_SOURCE_DIR}")
endif()

add_executable(ovm2ply ovm2ply/ovm2ply.cc)

target_link_libraries(ovm2ply
    OpenVolumeMesh::OpenVolumeMesh
    happly::happly
    )

set_target_properties (ovm2ply PROPERTIES
    CXX_STANDARD ${OVM_CXX_STANDARD}
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

#=]]
