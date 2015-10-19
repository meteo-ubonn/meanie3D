# -------------------------------------
# cmake configuration
# -------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_MACOSX_RPATH OFF)

# -------------------------------------
# Project configuration
# -------------------------------------

PROJECT(meanie3D-python)
GET_DIRECTORY_PROPERTY(PACKAGE_VERSION DIRECTORY .. DEFINITION PACKAGE_VERSION)

ADD_SUBDIRECTORY(meanie3D)