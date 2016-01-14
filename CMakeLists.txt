# CMAKE configuration file for project meanie3D
# (c) Jürgen Simon 2013 All Rights Reserved
# Email: juergen_simon@mac.com

# ------------------------------------- 
# cmake configuration
# -------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_MACOSX_RPATH OFF)

# ------------------------------------- 
# Project configuration
# -------------------------------------

PROJECT(meanie3D)
ENABLE_LANGUAGE(CXX)

SET(${PROJECT_NAME}_MAJOR_VERSION 1)
SET(${PROJECT_NAME}_MINOR_VERSION 6)
SET(${PROJECT_NAME}_PATCH_LEVEL 0)
SET(PACKAGE_VERSION "1.6.0")
SET(PROJECT_URL "http://git.meteo.uni-bonn.de/projects/meanie3d")
SET(PROJECT_LICENSE "MIT License")
SET(PROJECT_AUTHOR "Juergen Simon")
SET(PROJECT_AUTHOR_EMAIL "juergen.simon@uni-bonn.de")

# ------------------------------------- 
# Paths and build configuration
# -------------------------------------

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config_h.in
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}/config.h)

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/include/meanie3D/version.h.in
        ${PROJECT_SOURCE_DIR}/include/meanie3D/version.h)

# Set default cmake build type to Debug
# (None, Debug, Release, RelWithDebInfo or MinSizeRel)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF()

# Project include directory

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/include)

# ------------------------------------- 
# Options
# -------------------------------------

OPTION(WITH_OPENMP "Enable/Disable OpenMP" ON)
OPTION(WITH_TESTS "Enable/Disable building of tests" OFF)
OPTION(WITH_VTK "Enable/Disable visualisation code (requires VTK)" OFF)
OPTION(WITH_OPENCV "Legacy option. Leave at OFF" OFF)

# -------------------------------------
# Submodules
# -------------------------------------

ADD_SUBDIRECTORY(python)
ADD_SUBDIRECTORY(radolan)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/radolan/include)

# ------------------------------------- 
# Multithreading / Parallelism Choices
# -------------------------------------

SET(OPENMP_ENABLED ${WITH_OPENMP})
IF (OPENMP_ENABLED)
    MESSAGE(STATUS "OpenMP is enabled")
    ADD_DEFINITIONS(-DWITH_OPENMP=1)
ELSE()
    MESSAGE(STATUS "OpenMP is disabled")
ENDIF()

# -------------------------------------
# Visualization on/off
# -------------------------------------

SET(VTK_ENABLED ${WITH_VTK})
IF (VTK_ENABLED)
    MESSAGE(STATUS "VTK is enabled")
    ADD_DEFINITIONS(-DWITH_VTK=1)
ELSE()
    MESSAGE(STATUS "VTK is disabled")
ENDIF()

# -------------------------------------
# Automated tests yes/no
# -------------------------------------

SET(TESTS_ENABLED ${WITH_TESTS})
IF (TESTS_ENABLED)
    MESSAGE(STATUS "Tests are enabled.")
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(gtest-1.7.0)
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/gtest-1.7.0)
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/gtest-1.7.0/include)
ELSE()
    MESSAGE(STATUS "Tests are disabled. Use -DWITH_TESTS=ON/OFF to switch tests on or off.")
ENDIF()

# -------------------------------------
# OpenCV
# -------------------------------------

SET(OPENCV_ENABLED ${WITH_OPENCV})
IF (OPENCV_ENABLED)
    MESSAGE(STATUS "OpenCV is enabled.")
    ADD_DEFINITIONS(-DWITH_OPENMP=1)
ENDIF()

# -------------------------------------
# Compiler Flags
# -------------------------------------

MESSAGE(STATUS "CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")
MESSAGE(STATUS "CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "Compiler ID = ${CMAKE_CXX_COMPILER_ID}")
IF(CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE(STATUS "Compiler is GNU")
ENDIF()

IF (APPLE)
    IF (CMAKE_CXX_COMPILER MATCHES ".*clang")
        # Boost storage type problem
        ADD_DEFINITIONS(-DBOOST_HAS_INT128=1)
    ENDIF()
ENDIF()

IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pipe -ggdb -msse3 -g")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pipe -msse3")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -pipe -msse3")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -pipe -g -ggdb -msse3")
    IF(WITH_LTO)
        MESSAGE(STATUS "Link-Time Optimization switched on")
        SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
        SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -flto")
    ENDIF()
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ggdb -g")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ggdb")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -msse3 -O2")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -Os -msse3")
ENDIF()

# ------------------------------------- 
# Find dependencies
# -------------------------------------

MESSAGE(STATUS "Resolving dependencies for Meanie3D")

# netcdf
FIND_PACKAGE(NETCDF)
IF (NETCDF_FOUND)
    INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR})
ELSE(NETCDF_FOUND)
   MESSAGE(FATAL_ERROR "netcdf not found INCLUDE_DIR=${NETCDF_INCLUDE_DIR} LIBRARIES=${NETCDF_LIBRARIES}")
ENDIF(NETCDF_FOUND)

# HDF5
FIND_PACKAGE(HDF5)
IF (HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
ELSE(HDF5_FOUND)
   MESSAGE(FATAL_ERROR "HDF5 not found INCLUDE_DIR=${HDF5_INCLUDE_DIR} LIBRARIES=${HDF5_LIBRARIES}")
ENDIF(HDF5_FOUND)

# Boost
FIND_PACKAGE(Boost COMPONENTS program_options thread filesystem system)
IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ELSE(Boost_FOUND)
    MESSAGE(FATAL_ERROR ${Boost_ERROR_REASON})
ENDIF()

# FLANN
FIND_PACKAGE(FLANN)
IF (FLANN_FOUND)
    INCLUDE_DIRECTORIES(${FLANN_INCLUDE_DIR})
ELSE(FLANN_FOUND)
   MESSAGE(FATAL_ERROR "libflann not found (http://www.cs.ubc.ca/research/flann/)")
ENDIF(FLANN_FOUND)

# OpenMP
IF (OPENMP_ENABLED)
    FIND_PACKAGE(OpenMP REQUIRED)
    IF (OPENMP_FOUND)
        MESSAGE(STATUS "Found OpenMP (C_FLAGS=${OpenMP_C_FLAGS}, CXX_FLAGS=${OpenMP_CXX_FLAGS})")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE(OPENMP_FOUND)
        MESSAGE(ERROR "-- OpenMP is not enabled in your compiler")
    ENDIF (OPENMP_FOUND)

    IF (NOT CMAKE_COMPILER_IS_GNUCXX)
        FIND_PACKAGE(OpenMP_RT REQUIRED)
        IF (OpenMP_RT_FOUND)
            INCLUDE_DIRECTORIES(${OpenMP_RT_INCLUDE_DIRS})
        ELSE()
            MESSAGE(FATAL_ERROR "OpenMP Runtime not found (https://www.openmprtl.org)")
        ENDIF()
    ENDIF()

ENDIF(OPENMP_ENABLED)

# OpenCV
IF (OPENCV_ENABLED)
    FIND_PACKAGE(OpenCV)
    IF (OpenCV_FOUND)
        MESSAGE(STATUS "OpenCV FOUND")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenCV_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenCV_CXX_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenCV_EXE_LINKER_FLAGS}")
    ENDIF (OpenCV_FOUND)
ENDIF()

# VTK
IF (VTK_ENABLED)
    FIND_PACKAGE(VTK 6.0 REQUIRED NO_MODULE)
    IF (VTK_FOUND)
        MESSAGE(STATUS "VTK found")
        INCLUDE(${VTK_USE_FILE})
        INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
    ELSE()
        MESSAGE(FATAL_ERROR "VTK not found (http://www.vtk.org/VTK/resources/software.html)")
    ENDIF()
ENDIF()

# Blitz++
FIND_PACKAGE(Blitz REQUIRED)
IF (Blitz_FOUND)
    INCLUDE_DIRECTORIES(${Blitz_INCLUDE_DIR})
ELSE()
    MESSAGE(FATAL_ERROR "blitz++ not found")
ENDIF()

# ------------------------------------- 
# Binaries
# -------------------------------------

# Library

ADD_LIBRARY(meanie3D SHARED
    meanie3D/config.h	
    include/meanie3D/adaptors/Radolan.h
    include/meanie3D/adaptors.h
    include/meanie3D/array/array_index.h
    include/meanie3D/array/array_index_impl.h
    include/meanie3D/array/linear_index_mapping.h
    include/meanie3D/array/multiarray.h
    include/meanie3D/array/multiarray_blitz.h
    include/meanie3D/array/multiarray_boost.h
    include/meanie3D/array/multiarray_recursive.h
    include/meanie3D/array.h
    include/meanie3D/clustering/cluster.h
    include/meanie3D/clustering/cluster_impl.h
    include/meanie3D/clustering/cluster_list.h
    include/meanie3D/clustering/cluster_list_impl.h
    include/meanie3D/clustering/cluster_op.h
    include/meanie3D/clustering/cluster_op_impl.h
    include/meanie3D/clustering/cluster_utils.h
    include/meanie3D/clustering/cluster_utils_impl.h
    include/meanie3D/clustering/conrad_cluster.h
    include/meanie3D/clustering/histogram.h
    include/meanie3D/clustering/histogram_impl.h
    include/meanie3D/clustering/id.h
    include/meanie3D/clustering.h
    include/meanie3D/defines.h
    include/meanie3D/exceptions/CFFileConversionException.h
    include/meanie3D/exceptions.h
    include/meanie3D/featurespace/coordinate_system.h
    include/meanie3D/featurespace/coordinate_system_impl.h
    include/meanie3D/featurespace/data_store.h
    include/meanie3D/featurespace/featurespace.h
    include/meanie3D/featurespace/featurespace_impl.h
    include/meanie3D/featurespace/netcdf_data_store.h
    include/meanie3D/featurespace/point.h
    include/meanie3D/featurespace/point_default_factory.h
    include/meanie3D/featurespace/point_factory.h
    include/meanie3D/featurespace/point_impl.h
    include/meanie3D/featurespace/timestamp.h
    include/meanie3D/featurespace.h
    include/meanie3D/filters/convection_filter.h
    include/meanie3D/filters/convection_filter_impl.h
    include/meanie3D/filters/filter.h
    include/meanie3D/filters/scalespace_filter.h
    include/meanie3D/filters/scalespace_filter_impl.h
    include/meanie3D/filters/scalespace_kernel.h
    include/meanie3D/filters/scalespace_kernel_impl.h
    include/meanie3D/filters/replacement_filter.h
    include/meanie3D/filters/replacement_filter_impl.h
    include/meanie3D/filters/threshold_filter.h
    include/meanie3D/filters/threshold_filter_impl.h
    include/meanie3D/filters/weight_filter.h
    include/meanie3D/filters/weight_filter_impl.h
    include/meanie3D/filters.h
    include/meanie3D/implementations.h
    include/meanie3D/index/index.h
    include/meanie3D/index/index_impl.h
    include/meanie3D/index/index_whitening.h
    include/meanie3D/index/indexed_flann.h
    include/meanie3D/index/indexed_kdtree.h
    include/meanie3D/index/linear.h
    include/meanie3D/index/rectilinear_grid_index.h
    include/meanie3D/index/search_parameters.h
    include/meanie3D/index.h
    include/meanie3D/kdtree/kdtree.h
    include/meanie3D/meanie3D.h
    include/meanie3D/namespaces.h
    include/meanie3D/numericalrecipes/nrutil.h
    include/meanie3D/numericalrecipes/numericalrecipes.h
    include/meanie3D/operations/iterate_op.h
    include/meanie3D/operations/iterate_op_impl.h
    include/meanie3D/operations/kernels.h
    include/meanie3D/operations/kernels_impl.h
    include/meanie3D/operations/meanshift_op.h
    include/meanie3D/operations/meanshift_op_impl.h
    include/meanie3D/operations/operation.h
    include/meanie3D/operations.h
    include/meanie3D/parallel.h
    include/meanie3D/tracking/track.h
    include/meanie3D/tracking/track_cluster.h
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
    include/meanie3D/tracking.h
    include/meanie3D/utils/array_utils.h
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/file_utils.h
    include/meanie3D/utils/gaussian_normal.h
    include/meanie3D/utils/map_utils.h
    include/meanie3D/utils/matrix.h
    include/meanie3D/utils/matrix_impl.h
    include/meanie3D/utils/netcdf_utils.h
    include/meanie3D/utils/opencv_utils.h
    include/meanie3D/utils/rand_utils.h
    include/meanie3D/utils/set_utils.h
    include/meanie3D/utils/time_utils.h
    include/meanie3D/utils/vector_utils.h
    include/meanie3D/utils/verbosity.h
    include/meanie3D/utils/visit.h
    include/meanie3D/utils/visit_impl.h
    include/meanie3D/utils.h
    include/meanie3D/weights/brightband_evidence.h
    include/meanie3D/weights/ci_weights.h
    include/meanie3D/weights/default_weights.h
    include/meanie3D/weights/exp10_weight.h
    include/meanie3D/weights/inverse_default.h
    include/meanie3D/weights/oase_weights.h
    include/meanie3D/weights/weight_function.h
    include/meanie3D/weights.h
    src/adaptors/Radolan.cpp
    src/exceptions/CFFileConversionException.cpp
    src/kdtree/kdtree.c
    src/numericalrecipes/betai.c
    src/numericalrecipes/betacf.c
    src/numericalrecipes/crank.c
    src/numericalrecipes/erfcc.c
    src/numericalrecipes/gammln.c
    src/numericalrecipes/kendl1.c
    src/numericalrecipes/nrutil.c
    src/numericalrecipes/spear.c
    src/numericalrecipes/sort2.c
    src/utils/rand_utils.cpp
)

SOURCE_GROUP("meanie3d/adaptors" FILES
    include/meanie3D/adaptors/Radolan.h
)

SOURCE_GROUP("meanie3d/array" FILES
    include/meanie3D/array/array_index.h
    include/meanie3D/array/array_index_impl.h
    include/meanie3D/array/linear_index_mapping.h
    include/meanie3D/array/multiarray.h
    include/meanie3D/array/multiarray_blitz.h
    include/meanie3D/array/multiarray_boost.h
    include/meanie3D/array/multiarray_recursive.h
)

SOURCE_GROUP("meanie3d/clustering" FILES
    include/meanie3D/clustering/cluster.h
    include/meanie3D/clustering/cluster_impl.h
    include/meanie3D/clustering/cluster_list.h
    include/meanie3D/clustering/cluster_list_impl.h
    include/meanie3D/clustering/cluster_op.h
    include/meanie3D/clustering/cluster_op_impl.h
    include/meanie3D/clustering/cluster_utils.h
    include/meanie3D/clustering/cluster_utils_impl.h
    include/meanie3D/clustering/conrad_cluster.h
    include/meanie3D/clustering/histogram.h
    include/meanie3D/clustering/histogram_impl.h
    include/meanie3D/clustering/id.h
)

SOURCE_GROUP("meanie3d/exceptions" FILES
    include/meanie3D/exceptions/CFFileConversionException.h
)

SOURCE_GROUP("meanie3d/featurespace" FILES
    include/meanie3D/featurespace/coordinate_system.h
    include/meanie3D/featurespace/coordinate_system_impl.h
    include/meanie3D/featurespace/data_store.h
    include/meanie3D/featurespace/featurespace.h
    include/meanie3D/featurespace/featurespace_impl.h
    include/meanie3D/featurespace/netcdf_data_store.h
    include/meanie3D/featurespace/point.h
    include/meanie3D/featurespace/point_default_factory.h
    include/meanie3D/featurespace/point_factory.h
    include/meanie3D/featurespace/point_impl.h
    include/meanie3D/featurespace/timestamp.h
)

SOURCE_GROUP("meanie3d/filters" FILES
    include/meanie3D/filters/convection_filter.h
    include/meanie3D/filters/convection_filter_impl.h
    include/meanie3D/filters/filter.h
    include/meanie3D/filters/scalespace_filter.h
    include/meanie3D/filters/scalespace_filter_impl.h
    include/meanie3D/filters/scalespace_kernel.h
    include/meanie3D/filters/scalespace_kernel_impl.h
    include/meanie3D/filters/replacement_filter.h
    include/meanie3D/filters/replacement_filter_impl.h
    include/meanie3D/filters/threshold_filter.h
    include/meanie3D/filters/threshold_filter_impl.h
    include/meanie3D/filters/weight_filter.h
    include/meanie3D/filters/weight_filter_impl.h
)

SOURCE_GROUP("meanie3d/index" FILES
    include/meanie3D/index/index.h
    include/meanie3D/index/index_impl.h
    include/meanie3D/index/index_whitening.h
    include/meanie3D/index/indexed_flann.h
    include/meanie3D/index/indexed_kdtree.h
    include/meanie3D/index/linear.h
    include/meanie3D/index/rectilinear_grid_index.h
    include/meanie3D/index/search_parameters.h)

SOURCE_GROUP("meanie3d/operations" FILES
    include/meanie3D/operations/iterate_op.h
    include/meanie3D/operations/iterate_op_impl.h
    include/meanie3D/operations/kernels.h
    include/meanie3D/operations/kernels_impl.h
    include/meanie3D/operations/meanshift_op.h
    include/meanie3D/operations/meanshift_op_impl.h
    include/meanie3D/operations/operation.h
)

SOURCE_GROUP("meanie3d/tracking" FILES
    include/meanie3D/tracking/track.h
    include/meanie3D/tracking/track_cluster.h
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
)

SOURCE_GROUP("meanie3d/utils" FILES
    include/meanie3D/utils/array_utils.h
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/file_utils.h
    include/meanie3D/utils/gaussian_normal.h
    include/meanie3D/utils/map_utils.h
    include/meanie3D/utils/matrix.h
    include/meanie3D/utils/matrix_impl.h
    include/meanie3D/utils/netcdf_utils.h
    include/meanie3D/utils/opencv_utils.h
    include/meanie3D/utils/rand_utils.h
    include/meanie3D/utils/set_utils.h
    include/meanie3D/utils/time_utils.h
    include/meanie3D/utils/vector_utils.h
    include/meanie3D/utils/verbosity.h
    include/meanie3D/utils/visit.h
    include/meanie3D/utils/visit_impl.h
)

SOURCE_GROUP("meanie3d/weights" FILES
    include/meanie3D/weights/brightband_evidence.h
    include/meanie3D/weights/ci_weights.h
    include/meanie3D/weights/default_weights.h
    include/meanie3D/weights/exp10_weight.h
    include/meanie3D/weights/inverse_default.h
    include/meanie3D/weights/oase_weights.h
    include/meanie3D/weights/weight_function.h
)

# ------------------------------------------------------------- #
# Library and binaries
# ------------------------------------------------------------- #

TARGET_LINK_LIBRARIES(meanie3D
    radolan
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${Blitz_LIBRARY}
    ${OpenCV_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
#SET_TARGET_PROPERTIES(meanie3D PROPERTIES LINKER_LANGUAGE CXX)

# Conversion tool

ADD_EXECUTABLE(meanie3D-radolan2cfm
   src/executables/meanie3D-radolan2cfm.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-radolan2cfm
    meanie3D
    radolan
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-radolan2cfm PROPERTIES LINKER_LANGUAGE CXX)

# detection tool

ADD_EXECUTABLE(meanie3D-detect
    src/executables/meanie3D-detect.cpp)

TARGET_LINK_LIBRARIES(meanie3D-detect
    meanie3D
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-detect PROPERTIES LINKER_LANGUAGE CXX)

# tracking tool

ADD_EXECUTABLE(meanie3D-track
    meanie3D
    src/executables/meanie3D-track.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-track
    meanie3D
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-track PROPERTIES LINKER_LANGUAGE CXX)

# trackstats

ADD_EXECUTABLE(meanie3D-trackstats
    meanie3D
    src/executables/meanie3D-trackstats.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats
    meanie3D
    ${VTK_LIBRARIES}
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-trackstats PROPERTIES LINKER_LANGUAGE CXX)

# trackstats CONRAD

ADD_EXECUTABLE(meanie3D-trackstats-conrad
    meanie3D
    src/executables/meanie3D-trackstats-conrad.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats-conrad
    meanie3D
    ${VTK_LIBRARIES}
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-trackstats-conrad PROPERTIES LINKER_LANGUAGE CXX)

# copydims

ADD_EXECUTABLE(meanie3D-copydims
    meanie3D
    src/executables/meanie3D-copydims.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-copydims
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-copydims PROPERTIES LINKER_LANGUAGE CXX)

# Parallax correction for satellite 

ADD_EXECUTABLE(meanie3D-parallax_correction
    src/executables/meanie3D-parallax_correction.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-parallax_correction
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-parallax_correction PROPERTIES LINKER_LANGUAGE CXX)

# map data

ADD_EXECUTABLE(meanie3D-mapdata
    meanie3D
    src/executables/meanie3D-mapdata.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-mapdata
    meanie3D
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-mapdata PROPERTIES LINKER_LANGUAGE CXX)

# timestamp update

ADD_EXECUTABLE(meanie3D-timestamp
    src/executables/meanie3D-timestamp.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-timestamp
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-timestamp PROPERTIES LINKER_LANGUAGE CXX)

#
# VTK-converter
#
IF (VTK_ENABLED)
    ADD_EXECUTABLE(meanie3D-cfm2vtk
        src/executables/meanie3D-cfm2vtk.cpp)

    TARGET_LINK_LIBRARIES(meanie3D-cfm2vtk
        meanie3D
        ${Boost_LIBRARIES}
        ${VTK_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${OpenMP_RT_LIBRARIES})
    SET_TARGET_PROPERTIES(meanie3D-cfm2vtk PROPERTIES LINKER_LANGUAGE CXX)
ENDIF()

# valid_min/valid_max checker

ADD_EXECUTABLE(meanie3D-minmax
    src/executables/meanie3D-minmax.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-minmax
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-minmax PROPERTIES LINKER_LANGUAGE CXX)

ADD_EXECUTABLE(meanie3D-satconv
    src/executables/meanie3D-satconv.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-satconv
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-satconv PROPERTIES LINKER_LANGUAGE CXX)


# ------------------------------------------------------------- #
# Tests
# ------------------------------------------------------------- #

IF (TESTS_ENABLED)

    ADD_EXECUTABLE(m3D-test-kernel
        test/kernels/testcases.h
        test/kernels/test.cpp)

    TARGET_LINK_LIBRARIES(m3D-test-kernel
        meanie3D
        gtest
        ${Boost_LIBRARIES})
    SET_TARGET_PROPERTIES(m3D-test-kernel PROPERTIES LINKER_LANGUAGE CXX)

    # Unit tests for vector utils

    ADD_EXECUTABLE(m3D-test-collections
        test/collections/tests_arrayindex.h
        test/collections/tests_map.h
        test/collections/tests_multiarray.h
        test/collections/tests_set.h
        test/collections/tests_vector.h
        test/collections/test.cpp)

    TARGET_LINK_LIBRARIES(m3D-test-collections
        gtest
        meanie3D
        ${Boost_LIBRARIES})

    SET_TARGET_PROPERTIES(m3D-test-collections PROPERTIES LINKER_LANGUAGE CXX)

    # Unit tests for class FeatureSpace

    ADD_EXECUTABLE(m3D-test-featurespace
        test/settings.cpp
        test/settings.h
        test/testcase_base.h
        test/testcase_base_impl.h
        test/featurespace/circular_pattern.h
        test/featurespace/circular_pattern_impl.h
        test/featurespace/uniform.h
        test/featurespace/uniform_impl.h
        test/featurespace/weighed.h
        test/featurespace/weighed_impl.h
        test/featurespace/iteration.h
        test/featurespace/iteration_impl.h
        test/featurespace/testcases.h
        test/featurespace/test.cpp)

    TARGET_LINK_LIBRARIES(m3D-test-featurespace
        gtest
        meanie3D
        ${Boost_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${FLANN_LIBRARIES}
        ${OpenMP_RT_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES}
        ${Blitz_LIBRARY})

    SET_TARGET_PROPERTIES(m3D-test-featurespace PROPERTIES LINKER_LANGUAGE CXX)

    # KDTree

    ADD_EXECUTABLE(m3D-test-kdtree
        test/kdtree/test.c)

    TARGET_LINK_LIBRARIES(m3D-test-kdtree
        meanie3D
        ${Boost_LIBRARIES})
    SET_TARGET_PROPERTIES(m3D-test-kdtree PROPERTIES LINKER_LANGUAGE CXX)

    # Detection tests

    ADD_EXECUTABLE(m3D-test-detection
        test/settings.cpp
        test/settings.h
        test/testcase_base.h
        test/testcase_base_impl.h
        test/detection/clustering.h
        test/detection/clustering_impl.h
        test/detection/detection_test.h
        test/detection/test.cpp
        test/detection/testcases.h
        test/detection/variable_weighed.h
        test/detection/variable_weighed_impl.h)

    TARGET_LINK_LIBRARIES(m3D-test-detection
        gtest
        meanie3D
        ${Boost_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${OpenMP_RT_LIBRARIES}
        ${VTK_LIBRARIES})
    SET_TARGET_PROPERTIES(m3D-test-detection PROPERTIES LINKER_LANGUAGE CXX)

    ADD_TEST(KernelTest m3D-test-kernel)
    ADD_TEST(CollectionTest m3D-test-collections)
    ADD_TEST(KDTreeTest m3D-test-kdtree)
    ADD_TEST(MeanshiftTest m3D-test-featurespace)
    ADD_TEST(ClusteringTest m3D-test-detection)

    ADD_CUSTOM_TARGET(check COMMAND make test)

ENDIF(TESTS_ENABLED)

# ------------------------------------------------------------- #
# Installation
# ------------------------------------------------------------- #

INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/${PROJECT_NAME} 
    DESTINATION "/usr/local/include")

INSTALL(DIRECTORY include/${PROJECT_NAME}
        DESTINATION "/usr/local/include"
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN "*.swp" EXCLUDE)

INSTALL(TARGETS meanie3D LIBRARY DESTINATION "/usr/local/lib")
INSTALL(TARGETS meanie3D-detect RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-track RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-trackstats RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-trackstats-conrad RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-timestamp RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-parallax_correction RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-minmax RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-satconv RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-radolan2cfm RUNTIME DESTINATION "/usr/local/bin")
IF (WITH_VTK)
INSTALL(TARGETS meanie3D-cfm2vtk RUNTIME DESTINATION "/usr/local/bin")
ENDIF(WITH_VTK)

# ------------------------------------------------------------- #
# Documentation
# ------------------------------------------------------------- #

SET(FULL_VERSION_FOR_DOXYGEN 1.5.0)

INCLUDE(FindDoxygen)
IF(DOXYGEN_EXECUTABLE)

    IF(DOXYGEN_DOT_EXECUTABLE)
        SET(UBIQ_DOXYGEN_HAVE_DOT "YES")
    ELSE(DOXYGEN_DOT_EXECUTABLE)
        SET(UBIQ_DOXYGEN_HAVE_DOT "NO")
    ENDIF(DOXYGEN_DOT_EXECUTABLE)

    CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/doxy.in ${PROJECT_BINARY_DIR}/doxyfile.in)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/doc/html)
    ADD_CUSTOM_TARGET(docs  COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxyfile.in)

ENDIF(DOXYGEN_EXECUTABLE)
