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

# Python
FIND_PACKAGE(PythonInterp REQUIRED)
IF (PYTHONINTERP_FOUND)
    MESSAGE(STATUS "Found Python interpreter (${PYTHON_VERSION_STRING}) ${PYTHON_EXECUTABLE}")
ELSE()
    MESSAGE(FATAL_ERROR "python not found")
ENDIF()


# ------------------------------------------------------------- #
# Python packaging
# ------------------------------------------------------------- #
IF (PYTHONINTERP_FOUND)

    SET(SETUP_PY_IN "${CMAKE_CURRENT_BINARY_DIR}/setup.py.in")
    SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY})

    ADD_CUSTOM_TARGET(${PROJECT_NAME} ALL COMMAND ${PYTHON_EXECUTABLE} setup.py build)
    ADD_CUSTOM_TARGET(clean ALL COMMAND rm -rf build dist)

    INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${PYTHON_EXECUTABLE}\" ${SETUP_PY} install)")

ENDIF()
