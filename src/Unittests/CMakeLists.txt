SET(SOURCE_FILES
    unittests.cc
    unittests_attribs.cc
    unittests_basics.cc
    unittests_files.cc
    unittests_local_topology.cc
    unittests_common.cc
    unittests_handles.cc
    unittests_swap_entities.cc
    unittests_garbage_collection.cc
    unittests_iterators.cc
    unittests_mesh_copies.cc
    unittests_smart_tagger.cc
    unittests_properties.cc
    unittests_tet_mesh_navigation.cc
)

if(OVM_WITH_EIGEN)
    list(APPEND SOURCE_FILES unittests_eigen.cc)

    if (NOT TARGET Eigen3::Eigen AND OVM_DOWNLOAD_MISSING_EIGEN)
        include(FetchContent)
        set(FETCHCONTENT_QUIET OFF)
        set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)
        FetchContent_Declare(eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen
            GIT_TAG be56fff1ff288f83296a78d63314a3eeac3397ae # master 2025-10-26
            SOURCE_SUBDIR "nonexisting. Do not use Eigen CMake, we use it header-only."
        )
        FetchContent_MakeAvailable(eigen)
        message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")
        #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
    endif()
endif()

if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    find_package(OpenVolumeMesh REQUIRED)
endif()


include(FetchContent)

if (NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
    set (DEFINED FETCHCONTENT_UPDATES_DISCONNECTED ON)
endif()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.17.0
  EXCLUDE_FROM_ALL
)

if(WIN32)
    # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)



set(BUILD_GMOCK OFF)
include(GoogleTest)
set(INSTALL_GTEST OFF)

set_target_properties(gtest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

enable_testing()

set(TARGET_NAME "${OVM_TARGET_PREFIX}unittests")
add_executable(${TARGET_NAME} ${SOURCE_FILES})
target_link_libraries(${TARGET_NAME}
    OpenVolumeMesh::OpenVolumeMesh
    GTest::gtest GTest::gtest_main
    )

if(OVM_WITH_EIGEN)
    target_link_libraries(${TARGET_NAME} Eigen3::Eigen)
endif()

gtest_discover_tests(${TARGET_NAME}
    EXTRA_ARGS "--gtest_color=yes" "--gtest_output=xml"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests"
    )

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
    DESTINATION "${CMAKE_BINARY_DIR}/Unittests"
    )

#set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Unittests )

