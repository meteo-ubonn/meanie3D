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

ENABLE_LANGUAGE(CXX)

# ------------------------------------- 
# Project configuration
# -------------------------------------

PROJECT(meanie3D)

SET(${PROJECT_NAME}_MAJOR_VERSION 0)
SET(${PROJECT_NAME}_MINOR_VERSION 1)
SET(${PROJECT_NAME}_PATCH_LEVEL 0)

#setup OS

IF(CMAKE_SYSTEM_NAME MATCHES Darwin)
    SET(${PROJECT_NAME}_OS_MACOSX TRUE)
ELSE(CMAKE_SYSTEM_NAME MATCHES Darwin)
  IF(CMAKE_SYSTEM_NAME MATCHES Linux)
    SET(${PROJECT_NAME}_OS_LINUX TRUE)
  ELSE(CMAKE_SYSTEM_NAME MATCHES Linux)
    MESSAGE(FATAL_ERROR "unknown target OS")
  ENDIF(CMAKE_SYSTEM_NAME MATCHES Linux)
ENDIF(CMAKE_SYSTEM_NAME MATCHES Darwin)

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

# ------------------------------------- 
# Additional project configuration
# -------------------------------------

# Work around the problem in gtest with tuple
ADD_DEFINITIONS(-DGTEST_HAS_TR1_TUPLE=0 -DGTEST_HAVE_OWN_TR1_TUPLE=0)

# Switch on to use TBB
SET(WITH_TBB YES)

# Switch on to use Boost Threads. Do not combine with TBB
SET(WITH_BOOST_THREADS 0)

IF (WITH_TBB)
    ADD_DEFINITIONS(-DWITH_TBB=1)
ENDIF(WITH_TBB)

IF (WITH_BOOST_THREADS)
    ADD_DEFINITIONS(-DWITH_BOOST_THREADS=1)
ENDIF (WITH_BOOST_THREADS)

# -------------------------------------
# Compiler Flags
# -------------------------------------

SET(USE_CPP11 FALSE)

#Compiler flags

MESSAGE("CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")
MESSAGE("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
MESSAGE("Compiler ID = ${CMAKE_CXX_COMPILER_ID}")
IF(CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE("Compiler is G++")
ENDIF()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # Debug, Release, MinSizeRel and RelWithDebInfo
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -Wlocal-type-template-args -finline")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wlocal-type-template-args -flto -O3 -ftree-vectorize –vec-report2 -msse2 -finline")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -Wlocal-type-template-args -flto -Ofast -ftree-vectorize –vec-report2 -msse2 -finline")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wlocal-type-template-args -O1 -msse2 -finline")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    MESSAGE(ERROR "Compiler not currently supported")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    MESSAGE(ERROR "Compiler not currently supported")
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Debug, Release, MinSizeRel and RelWithDebInfo
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -finline -g -ggdb")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -flto")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -Os -flto")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O1 -g -ggdb")
ENDIF()

# C++11
IF (USE_CPP11)
   IF (APPLE)
      SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
      SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
   ELSE (APPLE)
      #TODO
   ENDIF (APPLE)
   SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++0x -stdlib=libc++ -g -Wall")
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
    ADD_DEFINITIONS(-DOPENMP=1)
    IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        MESSAGE(ERROR "OpenMP not supported for this compiler type")
    ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        MESSAGE(ERROR "OpenMP not supported for this compiler type")
    ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        MESSAGE(ERROR "OpenMP not supported for this compiler type")
    ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fopenmp")
        SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fopenmp")
        SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -fopenmp")
        SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fopenmp")
    ENDIF()
ENDIF (WITH_OPENMP)

# ------------------------------------- 
# Find dependencies
# -------------------------------------

MESSAGE(STATUS "Resolving dependencies for Meanie3D")

# netcdf
FIND_PACKAGE(NETCDF)
IF (NOT NETCDF_FOUND)
   MESSAGE(FATAL_ERROR "netcdf(cxx) not found (http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-cxx/)")
ENDIF()

# HDF5
FIND_PACKAGE(HDF5)
IF (NOT HDF5_FOUND)
   MESSAGE(FATAL_ERROR "HDF5 not found (http://www.hdfgroup.org/HDF5/)")
ENDIF()

# Boost
FIND_PACKAGE(Boost COMPONENTS system program_options thread filesystem)
IF(NOT Boost_FOUND)
    MESSAGE(FATAL_ERROR ${Boost_ERROR_REASON})
ENDIF()

# GTest
FIND_PACKAGE(GTest)
IF (NOT GTEST_FOUND)
   MESSAGE(FATAL_ERROR "libgtest not found (http://code.google.com/p/googletest/)")
ENDIF()

# TBB
IF (WITH_TBB)
    MESSAGE(STATUS "Looking for TBB")
    FIND_PACKAGE(TBB)
    IF (NOT TBB_FOUND)
        MESSAGE(FATAL_ERROR "TBB not found (http://threadingbuildingblocks.org)")
    ENDIF (NOT TBB_FOUND)
ENDIF(WITH_TBB)

IF (WITH_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
    IF (OPENMP_FOUND)
        MESSAGE("OPENMP FOUND")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF (OPENMP_FOUND)
ENDIF(WITH_OPENMP)

FIND_PACKAGE(OpenCV)
IF (OpenCV_FOUND)
    MESSAGE("OpenCV FOUND")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenCV_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenCV_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenCV_EXE_LINKER_FLAGS}")
ENDIF (OpenCV_FOUND)


# VTK
FIND_PACKAGE(VTK 6.0 REQUIRED NO_MODULE)
IF (NOT VTK_FOUND)
   MESSAGE(FATAL_ERROR "VTK not found (http://www.vtk.org/VTK/resources/software.html)")
ELSE(NOT VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ENDIF(NOT VTK_FOUND)

# Shapefile
FIND_PACKAGE(SHP)
IF (NOT SHP_FOUND)
   MESSAGE(FATAL_ERROR "Shapefile not found INCLUDE_DIR=${SHP_INCLUDE_DIR} LIBRARIES=${SHP_LIBRARIES}")
ENDIF()

# -------------------------------------
# Submodules
# -------------------------------------

ADD_SUBDIRECTORY(cf-algorithms)
ADD_SUBDIRECTORY(numericalrecipes)

# ------------------------------------- 
# Paths and Build Configuration
# -------------------------------------

IF(WITH_TBB)
	INCLUDE_DIRECTORIES(
		${PROJECT_SOURCE_DIR}/include
		${cf-algorithms_SOURCE_DIR}/include
		${cf-algorithms_SOURCE_DIR}/radolan/include
		${numericalrecipes_SOURCE_DIR}/include
		${NETCDF_INCLUDE_DIR}
		${HDF5_INCLUDE_DIR}
        ${SHP_INCLUDE_DIR}
        ${VTK_INCLUDE_DIRS}
        ${TBB_INCLUDE_DIRS}
        )
ELSE()
	INCLUDE_DIRECTORIES(
		${PROJECT_SOURCE_DIR}/include
		${cf-algorithms_SOURCE_DIR}/include
		${cf-algorithms_SOURCE_DIR}/radolan/include
		${numericalrecipes_SOURCE_DIR}/include
		${NETCDF_INCLUDE_DIR}
		${HDF5_INCLUDE_DIR}
        ${SHP_INCLUDE_DIR}
        ${VTK_INCLUDE_DIRS})
ENDIF()

CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/config_h.in 
	${PROJECT_BINARY_DIR}/${PROJECT_NAME}/config.h
)

#set default cmake build type to Debug

# (None Debug Release RelWithDebInfo MinSizeRel)
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF()

# ------------------------------------- 
# Binaries
# -------------------------------------

# Library

ADD_LIBRARY(meanie3D SHARED
    meanie3D/config.h	
    include/meanie3D/meanie3D.h
    include/meanie3D/defines.h
    include/meanie3D/namespaces.h
    include/meanie3D/types.h
    include/meanie3D/types/cluster.h
    include/meanie3D/types/cluster_impl.h
    include/meanie3D/types/cluster_list.h
    include/meanie3D/types/cluster_list_impl.h
    include/meanie3D/types/conrad_cluster.h
    include/meanie3D/types/point.h
    include/meanie3D/types/point_factory.h
    include/meanie3D/types/histogram.h
    include/meanie3D/types/histogram_impl.h
    include/meanie3D/types/timestamp.h
    include/meanie3D/filters.h
    include/meanie3D/filters/filter.h
    include/meanie3D/filters/scalespace_kernel.h
    include/meanie3D/filters/scalespace_kernel_impl.h
    include/meanie3D/filters/scalespace_filter.h
    include/meanie3D/filters/scalespace_filter_impl.h
    include/meanie3D/filters/threshold_filter.h
    include/meanie3D/filters/threshold_filter_impl.h
    include/meanie3D/filters/weight_filter.h
    include/meanie3D/filters/weight_filter_impl.h
    include/meanie3D/filters/convection_filter.h
    include/meanie3D/filters/convection_filter_impl.h
    include/meanie3D/detection.h
    include/meanie3D/detection/cluster_op.h
    include/meanie3D/detection/cluster_op_impl.h
    include/meanie3D/detection/cluster_task.h
    include/meanie3D/detection/cluster_utils.h
    include/meanie3D/detection/cluster_utils_impl.h
    include/meanie3D/tracking.h
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
    include/meanie3D/utils.h
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/verbosity.h
    include/meanie3D/utils/visit.h
    include/meanie3D/utils/visit_impl.h
    include/meanie3D/utils/matrix.h
    include/meanie3D/utils/matrix_impl.h
    include/meanie3D/utils/file_utils.h
    include/meanie3D/utils/opencv_utils.h
    include/meanie3D/weights/brightband_evidence.h
    include/meanie3D/weights/ci_weights.h
    include/meanie3D/weights/default_weights.h
    include/meanie3D/weights/exp10_weight.h
    include/meanie3D/weights/inverse_default.h
    include/meanie3D/weights/oase_weights.h
    include/meanie3D/weights/weights.h
    include/meanie3D/implementations.h
)

SOURCE_GROUP("meanie3d/types" FILES
    include/meanie3D/types/cluster.h
    include/meanie3D/types/cluster_impl.h
    include/meanie3D/types/cluster_list.h
    include/meanie3D/types/cluster_list_impl.h
    include/meanie3D/types/conrad_cluster.h
    include/meanie3D/types/point.h
    include/meanie3D/types/point_factory.h
    include/meanie3D/types/histogram.h
    include/meanie3D/types/histogram_impl.h
    include/meanie3D/types/timestamp.h
)

SOURCE_GROUP("meanie3d/filters" FILES
    include/meanie3D/filters/filter.h
    include/meanie3D/filters/scalespace_kernel.h
    include/meanie3D/filters/scalespace_kernel_impl.h
    include/meanie3D/filters/scalespace_filter.h
    include/meanie3D/filters/scalespace_filter_impl.h
    include/meanie3D/filters/threshold_filter.h
    include/meanie3D/filters/threshold_filter_impl.h
    include/meanie3D/filters/weight_filter.h
    include/meanie3D/filters/weight_filter_impl.h
    include/meanie3D/filters/convection_filter.h
    include/meanie3D/filters/convection_filter_impl.h
)

SOURCE_GROUP("meanie3d/detection" FILES
    include/meanie3D/detection/cluster_op.h
    include/meanie3D/detection/cluster_op_impl.h
    include/meanie3D/detection/cluster_task.h
    include/meanie3D/detection/cluster_utils.h
    include/meanie3D/detection/cluster_utils_impl.h
)

SOURCE_GROUP("meanie3d/tracking" FILES
    include/meanie3D/tracking/tracking.h
    include/meanie3D/tracking/tracking_impl.h
)

SOURCE_GROUP("meanie3d/utils" FILES
    include/meanie3D/utils/cluster_index.h
    include/meanie3D/utils/cluster_index_impl.h
    include/meanie3D/utils/file_utils.h
    include/meanie3D/utils/matrix.h
    include/meanie3D/utils/matrix_impl.h
    include/meanie3D/utils/opencv_utils.h
    include/meanie3D/utils/verbosity.h
    include/meanie3D/utils/visit.h
    include/meanie3D/utils/visit_impl.h
)

SOURCE_GROUP("meanie3d/weights" FILES
    include/meanie3D/weights/weights.h
    include/meanie3D/weights/ci_weights.h
    include/meanie3D/weights/oase_weights.h
    include/meanie3D/weights/exp10_weight.h
    include/meanie3D/weights/default_weights.h
    include/meanie3D/weights/inverse_default.h
    include/meanie3D/weights/brightband_evidence.h
)

TARGET_LINK_LIBRARIES(meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${BOOST_LIBRARIES}
    ${VTK_LIBRARIES}
)

# detection tool

ADD_EXECUTABLE(meanie3D-detect
    cf-algorithms
    numericalrecipes
    src/meanie3D-detect.cpp)

IF(WITH_TBB)
    TARGET_LINK_LIBRARIES(meanie3D-detect
        meanie3D
    	cf-algorithms
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${TBB_LIBRARIES}
        ${VTK_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3D-detect
        meanie3D
    	cf-algorithms
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES})
ENDIF()

# CI Score based detection

ADD_EXECUTABLE(meanie3D-ci-detect
    cf-algorithms
    numericalrecipes
    src/meanie3D-ci-detect.cpp)

IF(WITH_TBB)
    TARGET_LINK_LIBRARIES(meanie3D-ci-detect
        meanie3D
    	cf-algorithms
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${TBB_LIBRARIES}
        ${VTK_LIBRARIES}
        ${OpenCV_LIBS})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3D-ci-detect
        meanie3D
    	cf-algorithms
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES}
        ${OpenCV_LIBS})
ENDIF()

# tracking tool

ADD_EXECUTABLE(meanie3D-track
    cf-algorithms
    numericalrecipes
    src/meanie3D-track.cpp
)

IF(WITH_TBB)
    TARGET_LINK_LIBRARIES(meanie3D-track
        meanie3D
    	cf-algorithms
        numericalrecipes
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES}
        ${TBB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3D-track
        meanie3D
    	cf-algorithms
        numericalrecipes
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES})
ENDIF()

# trackstats

ADD_EXECUTABLE(meanie3D-trackstats
    cf-algorithms
    numericalrecipes
    src/meanie3D-trackstats.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})

# trackstats CONRAD

ADD_EXECUTABLE(meanie3D-trackstats-conrad
    cf-algorithms
    numericalrecipes
    src/meanie3D-trackstats-conrad.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-trackstats-conrad
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})

# copydims

ADD_EXECUTABLE(meanie3D-copydims
    cf-algorithms
    numericalrecipes
    src/meanie3D-copydims.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-copydims
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})
    
# Parallax correction for satellite 

ADD_EXECUTABLE(meanie3D-parallax_correction
    src/meanie3D-parallax_correction.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-parallax_correction
    meanie3D
    cf-algorithms
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES})

# map data

ADD_EXECUTABLE(meanie3D-mapdata
    cf-algorithms
    numericalrecipes
    src/meanie3D-mapdata.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-mapdata
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES}
    ${SHP_LIBRARIES})

# timestamp update

ADD_EXECUTABLE(meanie3D-timestamp
    cf-algorithms
    numericalrecipes
    src/meanie3D-timestamp.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-timestamp
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})

# VTK-converter

ADD_EXECUTABLE(meanie3D-cfm2vtk
    src/meanie3D-cfm2vtk.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-cfm2vtk
    meanie3D
    cf-algorithms
    ${Boost_LIBRARIES}
    ${VTK_LIBRARIES}
)

# valid_min/valid_max checker

ADD_EXECUTABLE(meanie3D-minmax
    src/meanie3D-minmax.cpp
)

TARGET_LINK_LIBRARIES(meanie3D-minmax
    meanie3D
    cf-algorithms
    numericalrecipes
    ${NETCDF_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${VTK_LIBRARIES})

# Unit Tests
# detection

ADD_EXECUTABLE(meanie3d-detection-test
#	src/test/gaussian_cloud.cpp
    src/test/detection/detection_test.h
    src/test/detection/settings.h
    src/test/detection/settings.cpp
    src/test/detection/testcase_base.h
    src/test/detection/testcase_base_impl.h
    src/test/detection/clustering.h
    src/test/detection/clustering_impl.h
    src/test/detection/testcases.h
    src/test/detection/test.cpp
)

IF(WITH_TBB)
    TARGET_LINK_LIBRARIES(meanie3d-detection-test
        meanie3D
    	cf-algorithms
        numericalrecipes
        ${Boost_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${GTEST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${VTK_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3d-detection-test
        meanie3D
    	cf-algorithms
        numericalrecipes
        ${Boost_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${VTK_LIBRARIES}
        ${GTEST_LIBRARIES})
ENDIF()

SET_TARGET_PROPERTIES(meanie3D PROPERTIES LINKER_LANGUAGE CXX)
SET_TARGET_PROPERTIES(meanie3D-detect PROPERTIES LINKER_LANGUAGE CXX)

# Installer 

INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION "/usr/local/include")
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

# Documentation

SET(FULL_VERSION_FOR_DOXYGEN 0.1.0)
