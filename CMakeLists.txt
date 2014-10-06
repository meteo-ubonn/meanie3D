# CMAKE configuration file for project meanie3D
# (c) Jürgen Simon 2013 All Rights Reserved
# Email: juergen_simon@mac.com

# ------------------------------------- 
# cmake configuration
# -------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_MACOSX_RPATH OFF)
ENABLE_LANGUAGE(CXX)

# ------------------------------------- 
# Project configuration
# -------------------------------------

PROJECT(meanie3D)

SET(${PROJECT_NAME}_MAJOR_VERSION 0)
SET(${PROJECT_NAME}_MINOR_VERSION 1)
SET(${PROJECT_NAME}_PATCH_LEVEL 0)

# ------------------------------------- 
# Paths and Build Configuration
# -------------------------------------

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/include)

CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/config_h.in 
	${PROJECT_BINARY_DIR}/${PROJECT_NAME}/config.h
)

#set default cmake build type to Debug

# (None Debug Release RelWithDebInfo MinSizeRel)
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF()

# ------------------------------------
#setup OS
# ------------------------------------

IF(CMAKE_SYSTEM_NAME MATCHES Darwin)
    SET(${PROJECT_NAME}_OS_MACOSX TRUE)
ELSE(CMAKE_SYSTEM_NAME MATCHES Darwin)
  IF(CMAKE_SYSTEM_NAME MATCHES Linux)
    SET(${PROJECT_NAME}_OS_LINUX TRUE)
  ELSE(CMAKE_SYSTEM_NAME MATCHES Linux)
    MESSAGE(FATAL_ERROR "unknown target OS")
  ENDIF(CMAKE_SYSTEM_NAME MATCHES Linux)
ENDIF(CMAKE_SYSTEM_NAME MATCHES Darwin)

# -------------------------------------
# Submodules
# -------------------------------------

ADD_SUBDIRECTORY(radolan)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/radolan/include)

# ------------------------------------- 
# Multithreading / Parallelism Choices
# -------------------------------------

IF (WITH_TBB)
    ADD_DEFINITIONS(-DWITH_TBB=1)
ENDIF(WITH_TBB)

IF (WITH_BOOST_THREADS)
    ADD_DEFINITIONS(-DWITH_BOOST_THREADS=1)
ENDIF (WITH_BOOST_THREADS)

IF (WITH_OPENMP)
    ADD_DEFINITIONS(-DWITH_OPENMP=1)
ENDIF()

# -------------------------------------
# Compiler Flags
# -------------------------------------

MESSAGE("-- CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")
MESSAGE("-- CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
MESSAGE("-- Compiler ID = ${CMAKE_CXX_COMPILER_ID}")
IF(CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE("-- Compiler is G++")
ENDIF()

IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -finline")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ftree-vectorize -msse2 -finline")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -ftree-vectorize -msse2 -finline")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ftree-vectorize -msse2 -finline")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ggdb")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ggdb")
ENDIF()

IF(WITH_LTO)
    MESSAGE(STATUS "Link-Time Optimization switched on")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -flto")
ELSE()
    
ENDIF()

# C++11
IF (USE_CPP11)
   IF (APPLE)
      SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
      SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
   ELSE (APPLE)
      #TODO
   ENDIF (APPLE)
   SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++0x -stdlib=libc++ -g")
   SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")
   SET(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++")
ELSE(USE_CPP11)
   IF (APPLE)
      SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wc++11-extensions ")
   ELSE (APPLE)
   ENDIF (APPLE)
ENDIF (USE_CPP11)

# OpenMP 

IF (WITH_OPENMP)
    IF (CMAKE_CXX_COMPILER MATCHES ".*clang")
        ADD_DEFINITIONS(-DBOOST_HAS_INT128=1)
    ENDIF()
    ADD_DEFINITIONS(-DOPENMP=1)
ENDIF (WITH_OPENMP)

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

# GTest
FIND_PACKAGE(GTest)
IF(GTEST_FOUND)
    INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIR})
    ADD_DEFINITIONS(-DGTEST_HAS_TR1_TUPLE=0 -DGTEST_HAVE_OWN_TR1_TUPLE=0)
ELSE(GTEST_FOUND)
    MESSAGE(FATAL_ERROR "libgtest not found INCLUDE_DIR=${GTEST_INCLUDE_DIR} LIBRARIES=${GTEST_LIBRARIES}")
    MESSAGE(FATAL_ERROR "did you specify -D GTEST_INCLUDE_DIR and -D GTEST_LIBRARIES ?")
ENDIF(GTEST_FOUND)

# FLANN
FIND_PACKAGE(FLANN)
IF (FLANN_FOUND)
    INCLUDE_DIRECTORIES(${FLANN_INCLUDE_DIR})
ELSE(FLANN_FOUND)
   MESSAGE(FATAL_ERROR "libflann not found (http://www.cs.ubc.ca/research/flann/)")
ENDIF(FLANN_FOUND)

# TBB
IF (WITH_TBB)
    FIND_PACKAGE(TBB)
    IF (TBB_FOUND)
        INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
    ELSE(TBB_FOUND)
        MESSAGE(FATAL_ERROR "TBB not found (http://threadingbuildingblocks.org)")
   ENDIF(TBB_FOUND)
ENDIF(WITH_TBB)

# OpenMP
IF (WITH_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
    IF (OPENMP_FOUND)
        MESSAGE(STATUS "Found OpenMP (C_FLAGS=${OpenMP_C_FLAGS}, CXX_FLAGS=${OpenMP_CXX_FLAGS})")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE(OPENMP_FOUND)
        MESSAGE(ERROR "-- OpenMP is not enabled in your compiler")
    ENDIF (OPENMP_FOUND)

    FIND_PACKAGE(OpenMP_RT REQUIRED)
    IF (OpenMP_RT_FOUND)
        INCLUDE_DIRECTORIES(${OpenMP_RT_INCLUDE_DIRS})
    ELSE(OpenMP_RT_FOUND)
        MESSAGE(FATAL_ERROR "OpenMP Runtime not found (https://www.openmprtl.org)")
    ENDIF(OpenMP_RT_FOUND)

ENDIF(WITH_OPENMP)

# OpenCV
FIND_PACKAGE(OpenCV)
IF (OpenCV_FOUND)
    MESSAGE("-- OpenCV FOUND")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenCV_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenCV_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenCV_EXE_LINKER_FLAGS}")
ENDIF (OpenCV_FOUND)

# VTK
FIND_PACKAGE(VTK 6.0 REQUIRED NO_MODULE)
IF (VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ELSE()
    MESSAGE(FATAL_ERROR "VTK not found (http://www.vtk.org/VTK/resources/software.html)")
ENDIF()

# Shapefile
FIND_PACKAGE(SHP)
IF (SHP_FOUND)
    INCLUDE_DIRECTORIES(${SHP_INCLUDE_DIRS})
ELSE()
   MESSAGE(FATAL_ERROR "Shapefile not found INCLUDE_DIR=${SHP_INCLUDE_DIR} LIBRARIES=${SHP_LIBRARIES}")
ENDIF()

# Blitz++
FIND_PACKAGE(Blitz REQUIRED)
IF (Blitz_FOUND)
    INCLUDE_DIRECTORIES(${Blitz_INCLUDE_DIR})
ELSE(Blitz_FOUND)
    MESSAGE(FATAL_ERROR "blitz++ not found")
ENDIF(Blitz_FOUND)

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
    include/meanie3D/clustering/cluster_task.h
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
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
    include/meanie3D/tracking.h
    include/meanie3D/utils/array_utils.h
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/console_utils.h
    include/meanie3D/utils/console_utils_impl.h
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
    include/meanie3D/clustering/cluster_task.h
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
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
)

SOURCE_GROUP("meanie3d/utils" FILES
    include/meanie3D/utils/array_utils.h
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/console_utils.h
    include/meanie3D/utils/console_utils_impl.h
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
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
    ${TBB_LIBRARIES}
    ${OpenMP_RT_LIBRARIES}
    ${Blitz_LIBRARY}
    ${OpenCV_LIBRARIES}
)
SET_TARGET_PROPERTIES(meanie3D PROPERTIES LINKER_LANGUAGE CXX)

# Conversion tool

ADD_EXECUTABLE(meanie3D-radolan2cfm
   src/executables/meanie3D-radolan2cfm.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-radolan2cfm
    radolan
    meanie3D
    ${RADOLAN_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-radolan2cfm PROPERTIES LINKER_LANGUAGE CXX)

# detection tool

ADD_EXECUTABLE(meanie3D-detect
    src/executables/meanie3D-detect.cpp)

TARGET_LINK_LIBRARIES(meanie3D-detect
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${TBB_LIBRARIES}
    ${VTK_LIBRARIES}
    ${OpenMP_RT_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-detect PROPERTIES LINKER_LANGUAGE CXX)

# tracking tool

ADD_EXECUTABLE(meanie3D-track
    meanie3D
    src/executables/meanie3D-track.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-track
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES}
    ${OpenMP_RT_LIBRARIES}
    ${TBB_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-track PROPERTIES LINKER_LANGUAGE CXX)

# trackstats

ADD_EXECUTABLE(meanie3D-trackstats
    meanie3D
    src/executables/meanie3D-trackstats.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-trackstats PROPERTIES LINKER_LANGUAGE CXX)

# trackstats CONRAD

ADD_EXECUTABLE(meanie3D-trackstats-conrad
    meanie3D
    src/executables/meanie3D-trackstats-conrad.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats-conrad
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-trackstats-conrad PROPERTIES LINKER_LANGUAGE CXX)

# copydims

ADD_EXECUTABLE(meanie3D-copydims
    meanie3D
    src/executables/meanie3D-copydims.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-copydims
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-copydims PROPERTIES LINKER_LANGUAGE CXX)

# Parallax correction for satellite 

ADD_EXECUTABLE(meanie3D-parallax_correction
    src/executables/meanie3D-parallax_correction.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-parallax_correction
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-parallax_correction PROPERTIES LINKER_LANGUAGE CXX)

# map data

ADD_EXECUTABLE(meanie3D-mapdata
    meanie3D
    src/executables/meanie3D-mapdata.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-mapdata
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES}
    ${SHP_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-mapdata PROPERTIES LINKER_LANGUAGE CXX)

# timestamp update

ADD_EXECUTABLE(meanie3D-timestamp
    src/executables/meanie3D-timestamp.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-timestamp
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-timestamp PROPERTIES LINKER_LANGUAGE CXX)

# VTK-converter

ADD_EXECUTABLE(meanie3D-cfm2vtk
    src/executables/meanie3D-cfm2vtk.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-cfm2vtk
    meanie3D
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-cfm2vtk PROPERTIES LINKER_LANGUAGE CXX)

# valid_min/valid_max checker

ADD_EXECUTABLE(meanie3D-minmax
    src/executables/meanie3D-minmax.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-minmax
    meanie3D
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(meanie3D-minmax PROPERTIES LINKER_LANGUAGE CXX)

# ------------------------------------------------------------- #
# Tests
# ------------------------------------------------------------- #

ADD_EXECUTABLE(m3D-test-kernel
    test/kernels/testcases.h
    test/kernels/test.cpp)

TARGET_LINK_LIBRARIES(m3D-test-kernel
    meanie3D
    ${Boost_LIBRARIES}
    ${GTEST_LIBRARIES})
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
    meanie3D
    ${Boost_LIBRARIES}
    ${GTEST_LIBRARIES})
SET_TARGET_PROPERTIES(m3D-test-collections PROPERTIES LINKER_LANGUAGE CXX)

# Unit tests for class FeatureSpace

ADD_EXECUTABLE(m3D-test-featurespace
    test/featurespace/settings.h
    test/featurespace/settings.cpp
    test/featurespace/testcase_base.h
    test/featurespace/testcase_base_impl.h
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
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${FLANN_LIBRARIES}
    ${TBB_LIBRARIES}
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
    test/detection/clustering.h
    test/detection/clustering_impl.h
    test/detection/detection_test.h
    test/detection/settings.cpp
    test/detection/settings.h
    test/detection/test.cpp
    test/detection/testcase_base.h
    test/detection/testcase_base_impl.h
    test/detection/testcases.h
    test/detection/variable_weighed.h
    test/detection/variable_weighed_impl.h)

TARGET_LINK_LIBRARIES(m3D-test-detection
    meanie3D
    ${Boost_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${TBB_LIBRARIES}
    ${OpenMP_RT_LIBRARIES}
    ${VTK_LIBRARIES})
SET_TARGET_PROPERTIES(m3D-test-detection PROPERTIES LINKER_LANGUAGE CXX)

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
INSTALL(TARGETS meanie3D-cfm2vtk RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-parallax_correction RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-minmax RUNTIME DESTINATION "/usr/local/bin")
INSTALL(TARGETS meanie3D-radolan2cfm RUNTIME DESTINATION "/usr/local/bin")

# ------------------------------------------------------------- #
# Documentation
# ------------------------------------------------------------- #

SET(FULL_VERSION_FOR_DOXYGEN 0.1.0)

INCLUDE(FindDoxygen)
IF(DOXYGEN_EXECUTABLE)

    IF(DOXYGEN_DOT_EXECUTABLE)
        SET(UBIQ_DOXYGEN_HAVE_DOT "YES")
    ELSE(DOXYGEN_DOT_EXECUTABLE)
        SET(UBIQ_DOXYGEN_HAVE_DOT "NO")
    ENDIF(DOXYGEN_DOT_EXECUTABLE)

    CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/doxy.in
                   ${PROJECT_BINARY_DIR}/doxyfile)

    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/doc/html)

    ADD_CUSTOM_TARGET(documentation 
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxyfile)

ENDIF(DOXYGEN_EXECUTABLE)
