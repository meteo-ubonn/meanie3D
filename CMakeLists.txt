# CMAKE configuration file for project meanie3D
# (c) Jürgen Simon 2013 All Rights Reserved
# Email: juergen_simon@mac.com

# ------------------------------------- 
# cmake configuration
# -------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
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
SET(WITH_TBB NO)

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

SET(USE_GCC FALSE)
SET(USE_CPP11 FALSE)

#enable auto-vectorization
IF (USE_GCC)
    SET(CMAKE_CXX_COMPILER=$CXX)
    IF (APPLE)
        SET(CMAKE_XCODE_ATTRIBUTE_GCC_VERSION "com.apple.compilers.llvmgcc42")
    ELSE (APPLE)
      #TODO
    ENDIF (APPLE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -msse2 -ffast-math")
ELSE (USE_GCC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -mllvm -vectorize -ffast-math")
ENDIF (USE_GCC)

# Optimization flags
IF (APPLE)
   SET(CMAKE_XCODE_ATTRIBUTE_LLVM_LTO[variant=Release] YES)
   SET(CMAKE_XCODE_ATTRIBUTE_GCC_STRICT_ALIASING[variant=Release] YES)
   SET(CMAKE_XCODE_ATTRIBUTE_GCC_UNROLL_LOOPS[variant=Release] YES)
   SET(CMAKE_XCODE_ATTRIBUTE_GCC_FAST_MATH[variant=Release] YES)   
ELSE (APPLE)
   #TODO
ENDIF(APPLE)   

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
ENDIF (USE_CPP11)

# ------------------------------------- 
# Find dependencies
# -------------------------------------

# netcdf
FIND_PACKAGE(NETCDF)
IF (NOT NETCDF_FOUND)
   MESSAGE(FATAL_ERROR "netcdf not found INCLUDE_DIR=${NETCDF_INCLUDE_DIR} LIBRARIES=${NETCDF_LIBRARIES}")
ENDIF()

# HDF5
FIND_PACKAGE(HDF5)
IF (NOT HDF5_FOUND)
   MESSAGE(FATAL_ERROR "HDF5 not found INCLUDE_DIR=${HDF5_INCLUDE_DIR} LIBRARIES=${HDF5_LIBRARIES}")
ENDIF()

# Boost
FIND_PACKAGE(Boost COMPONENTS system program_options thread filesystem)
IF(NOT Boost_FOUND)
    MESSAGE(FATAL_ERROR ${Boost_ERROR_REASON})
ENDIF()

# GTest
FIND_PACKAGE(GTest)
IF (NOT GTEST_FOUND)
   MESSAGE(FATAL_ERROR "libgtest not found INCLUDE_DIR=${GTEST_INCLUDE_DIR} LIBRARIES=${GTEST_LIBRARIES}")
ENDIF()

# TBB
IF (WITH_TBB)
   FIND_PACKAGE(TBB)
   IF (NOT TBB_FOUND)
      MESSAGE(FATAL_ERROR "TBB not found (http://threadingbuildingblocks.org)")
   ENDIF(NOT TBB_FOUND)
ENDIF(WITH_TBB)

# cf-algorithms
FIND_PACKAGE(CFA)
IF (NOT CFA_FOUND)
   MESSAGE(FATAL_ERROR "cf-algorithms not found INCLUDE_DIR=${CFA_INCLUDE_DIR} LIBRARIES=${CFA_LIBRARIES}")
ENDIF()

# ------------------------------------- 
# Paths and Build Configuration
# -------------------------------------

IF(WITH_TBB)
	INCLUDE_DIRECTORIES(
		${PROJECT_SOURCE_DIR}/include
		${NETCDF_INCLUDE_DIR}
		${HDF5_INCLUDE_DIR}
		${CFA_INCLUDE_DIR}
		${RADOLAN_INCLUDE_DIR}
		${TBB_INCLUDE_DIRS})
ELSE()
	INCLUDE_DIRECTORIES(
		${PROJECT_SOURCE_DIR}/include
		${NETCDF_INCLUDE_DIR}
		${HDF5_INCLUDE_DIR}
		${CFA_INCLUDE_DIR}
		${RADOLAN_INCLUDE_DIR})
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
    include/meanie3D/types/point.h
    include/meanie3D/types/point_impl.h
	include/meanie3D/filters.h
    include/meanie3D/filters/filter.h
    include/meanie3D/filters/scalespace_filter.h
    include/meanie3D/filters/scalespace_filter_impl.h
    include/meanie3D/filters/threshold_filter.h
    include/meanie3D/filters/threshold_filter_impl.h
	include/meanie3D/detection.h
    include/meanie3D/detection/cluster_op.h
    include/meanie3D/detection/cluster_op_impl.h
    include/meanie3D/detection/cluster_task.h
    include/meanie3D/utils.h
    include/meanie3D/utils/visit.h
	include/meanie3D/implementations.h
)

TARGET_LINK_LIBRARIES(meanie3D
   ${NETCDF_LIBRARIES}
   ${CFA_LIBRARIES}
)

# Executables

# clustering tool
ADD_EXECUTABLE(meanie3D-detect
    src/meanie3D-detect.cpp
)

IF(WITH_TBB)
    TARGET_LINK_LIBRARIES(meanie3D-detect
        meanie3D
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${CFA_LIBRARIES}
		${TBB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3D-detect
        meanie3D
        ${NETCDF_LIBRARIES}
        ${Boost_LIBRARIES}
        ${HDF5_LIBRARIES}
		${CFA_LIBRARIES})
ENDIF()

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
        ${Boost_LIBRARIES}
        ${NETCDF_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${GTEST_LIBRARIES}
		${CFA_LIBRARIES}
        ${TBB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(meanie3d-detection-test
	    meanie3D
	    ${Boost_LIBRARIES}
	    ${NETCDF_LIBRARIES}
	    ${HDF5_LIBRARIES}
	    ${GTEST_LIBRARIES}
		${CFA_LIBRARIES}
	    ${TBB_LIBRARIES})
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

# Documentation

SET(FULL_VERSION_FOR_DOXYGEN 0.1.0)