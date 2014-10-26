# - Find python
# This module searches for both the python interpreter and the python libraries
# and determines where they are located
#
#
#  PYTHON_FOUND - The requested Python components were found
#  PYTHON_EXECUTABLE  - path to the Python interpreter
#  PYTHON_INCLUDE_DIRS - path to the Python header files
#  PYTHON_LIBRARIES - the libraries to link against for python
#

#=============================================================================
# Copyright 2010 Branan Purvine-Riley
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

IF("3" STREQUAL "${Python_FIND_VERSION_MAJOR}")
  SET(PYTHON_3_OK "TRUE")
  SET(PYTHON_2_OK "FALSE") # redundant in version selection code, but skips a FOREACH
ELSE("3" STREQUAL "${Python_FIND_VERSION_MAJOR}")
  SET(PYTHON_2_OK "TRUE")
  # don't set PYTHON_3_OK to false here - if the user specified it we want to search for Python 2 & 3
ENDIF("3" STREQUAL "${Python_FIND_VERSION_MAJOR}")

# This is  heavily inspired by FindBoost.cmake, with the addition of a second version list to keep
# python 2 and python 3 separate
IF(Python_FIND_VERSION_EXACT)
  SET(_PYTHON_TEST_VERSIONS "${Python_FIND_VERSION_MAJOR}.${Python_FIND_VERSION_MINOR}")
ELSE(Python_FIND_VERSION_EXACT)
  SET(_PYTHON_3_KNOWN_VERSIONS ${PYTHON_3_ADDITIONAL_VERSIONS}
    "3.1" "3.0")
  SET(_PYTHON_2_KNOWN_VERSIONS ${PYTHON_2_ADDITIONAL_VERSIONS}
    "2.7" "2.6" "2.5" "2.4" "2.3" "2.2" "2.1" "2.0" "1.6" "1.5")
  SET(_PYTHON_TEST_VERSIONS)
  IF(Python_FIND_VERSION)
    IF(PYTHON_3_OK)
      FOREACH(version ${_PYTHON_3_KNOWN_VERSIONS})
        IF(NOT ${version} VERSION_LESS ${Python_FIND_VERSION})
          LIST(APPEND _PYTHON_TEST_VERSIONS ${version})
        ENDIF(NOT ${version} VERSION_LESS ${Python_FIND_VERSION})
      ENDFOREACH(version)
    ENDIF(PYTHON_3_OK)
    IF(PYTHON_2_OK)
      FOREACH(version ${_PYTHON_2_KNOWN_VERSIONS})
        IF(NOT ${version} VERSION_LESS ${Python_FIND_VERSION})
          LIST(APPEND _PYTHON_TEST_VERSIONS ${version})
        ENDIF(NOT ${version} VERSION_LESS ${Python_FIND_VERSION})
      ENDFOREACH(version)
    ENDIF(PYTHON_2_OK)
  ELSE(Python_FIND_VERSION)
    IF(PYTHON_3_OK)
      LIST(APPEND _PYTHON_TEST_VERSIONS ${_PYTHON_3_KNOWN_VERSIONS})
    ENDIF(PYTHON_3_OK)
    IF(PYTHON_2_OK)
      LIST(APPEND _PYTHON_TEST_VERSIONS ${_PYTHON_2_KNOWN_VERSIONS})
    ENDIF(PYTHON_2_OK)
  ENDIF(Python_FIND_VERSION)
ENDIF(Python_FIND_VERSION_EXACT)

SET(_PYTHON_EXE_VERSIONS)
FOREACH(version ${_PYTHON_TEST_VERSIONS})
  LIST(APPEND _PYTHON_EXE_VERSIONS python${version})
ENDFOREACH(version ${_PYTHON_TEST_VERSIONS})

IF(WIN32)
  SET(_PYTHON_REGISTRY_KEYS)
  FOREACH(version ${_PYTHON_TEST_VERSIONS})
    LIST(APPEND _PYTHON_REGISTRY_KEYS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\${version}\\InstallPath])
  ENDFOREACH(version ${_PYTHON_TEST_VERSIONS})
  # this will find any standard windows Python install before it finds anything from Cygwin
  FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python ${_PYTHON_EXE_VERSIONS} PATHS ${_PYTHON_REGISTRY_KEYS})
ELSE(WIN32)
  FIND_PROGRAM(PYTHON_EXECUTABLE NAMES ${_PYTHON_EXE_VERSIONS} python PATHS)
ENDIF(WIN32)

EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "-c" "from sys import *; stdout.write(str(version_info[0])+'.'+str(version_info[1]))" OUTPUT_VARIABLE _PYTHON_VERSION)


# Make sure our python version matches the requested version when we get a non-versioned executable name
GET_FILENAME_COMPONENT(_PYTHON_EXENAME ${PYTHON_EXECUTABLE} NAME_WE)
IF(Python_FIND_VERSION AND PYTHON_EXECUTABLE AND ${_PYTHON_EXENAME} STREQUAL "python")
  STRING(SUBSTRING "${_PYTHON_VERSION}" 7 3 _PYTHON_VERSION)
  IF(Python_FIND_VERSION_EXACT)
    IF(NOT ${_PYTHON_VERSION} VERSION_EQUAL "${Python_FIND_VERSION_MAJOR}.${Python_FIND_VERSION_MINOR}")
      SET(PYTHON_EXECUTABLE PYTHON_EXECUTABLE-NOTFOUND)
    ENDIF(NOT ${_PYTHON_VERSION} VERSION_EQUAL "${Python_FIND_VERSION_MAJOR}.${Python_FIND_VERSION_MINOR}")
  ELSE(Python_FIND_VERSION_EXACT)
    IF(NOT ${_PYTHON_VERSION} VERSION_LESS ${Python_FIND_VERSION})
      SET(PYTHON_EXECUTABLE PYTHON_EXECUTABLE-NOTFOUND)
    ENDIF(NOT ${_PYTHON_VERSION} VERSION_LESS ${Python_FIND_VERSION})
  ENDIF(Python_FIND_VERSION_EXACT)
ENDIF(Python_FIND_VERSION AND PYTHON_EXECUTABLE AND ${_PYTHON_EXENAME} STREQUAL "python")

# Get our library path and include directory from python
IF(WIN32)
  EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" # maybe there's a simpler way to do this? I haven't found it yet
    "-c" "import distutils.sysconfig; import sys; import os; sys.stdout.write(distutils.sysconfig.get_config_var('LIBRARY'))"
    OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
  )
ELSE(WIN32)
  EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" # maybe there's a simpler way to do this? I haven't found it yet
    "-c" "import distutils.sysconfig; import sys; import os; sys.stdout.write(distutils.sysconfig.get_config_var('LDLIBRARY'))"
    OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
  )
ENDIF(WIN32)
EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" # maybe there's a simpler way to do this? I haven't found it yet
  "-c" "import distutils.sysconfig; import sys; import os; sys.stdout.write(distutils.sysconfig.get_config_var('LIBDIR'))"
  OUTPUT_VARIABLE _PYTHON_LIBRARY_PATH
  )
EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" 
    "-c" "import distutils.sysconfig; import sys; import os; sys.stdout.write(distutils.sysconfig.get_config_var('MULTIARCH') or '')"
  OUTPUT_VARIABLE _PYTHON_LIBRARY_MULTIARCH
)
EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}"
  "-c" "import distutils.sysconfig; import sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
  OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
)

FIND_FILE(PYTHON_LIBRARY 
    ${_PYTHON_LIBRARY_NAME} 
    PATHS 
        ${_PYTHON_LIBRARY_PATH} 
        ${_PYTHON_LIBRARY_PATH}/${_PYTHON_LIBRARY_MULTIARCH}
    NO_DEFAULT_PATH
    )
FIND_FILE(PYTHON_HEADER "Python.h" PATHS ${PYTHON_INCLUDE_DIR} NO_DEFAULT_PATH)


set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR})
set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Python DEFAULT_MSG PYTHON_EXECUTABLE PYTHON_HEADER PYTHON_LIBRARY)

set (PYTHON_VERSION ${_PYTHON_VERSION})
MARK_AS_ADVANCED(PYTHON_EXECUTABLE)
MARK_AS_ADVANCED(PYTHON_INCLUDE_DIRS)
MARK_AS_ADVANCED(PYTHON_LIBRARIES)
MARK_AS_ADVANCED(PYTHON_INCLUDE_DIR)
MARK_AS_ADVANCED(PYTHON_LIBRARY)
MARK_AS_ADVANCED(PYTHON_VERSION)



# PYTHON_ADD_MODULE(<name> src1 src2 ... srcN) is used to build modules for python.
# PYTHON_WRITE_MODULES_HEADER(<filename>) writes a header file you can include
# in your sources to initialize the static python modules
FUNCTION(PYTHON_ADD_MODULE _NAME )
  GET_PROPERTY(_TARGET_SUPPORTS_SHARED_LIBS
    GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
  OPTION(PYTHON_ENABLE_MODULE_${_NAME} "Add module ${_NAME}" TRUE)
  OPTION(PYTHON_MODULE_${_NAME}_BUILD_SHARED
    "Add module ${_NAME} shared" ${_TARGET_SUPPORTS_SHARED_LIBS})

  # Mark these options as advanced
  MARK_AS_ADVANCED(PYTHON_ENABLE_MODULE_${_NAME}
    PYTHON_MODULE_${_NAME}_BUILD_SHARED)

  IF(PYTHON_ENABLE_MODULE_${_NAME})
    IF(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      SET(PY_MODULE_TYPE MODULE)
    ELSE(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      SET(PY_MODULE_TYPE STATIC)
      SET_PROPERTY(GLOBAL  APPEND  PROPERTY  PY_STATIC_MODULES_LIST ${_NAME})
    ENDIF(PYTHON_MODULE_${_NAME}_BUILD_SHARED)

    SET_PROPERTY(GLOBAL  APPEND  PROPERTY  PY_MODULES_LIST ${_NAME})
    ADD_LIBRARY(${_NAME} ${PY_MODULE_TYPE} ${ARGN})
#    TARGET_LINK_LIBRARIES(${_NAME} ${PYTHON_LIBRARIES})

    IF(PYTHON_MODULE_${_NAME}_BUILD_SHARED)
      SET_TARGET_PROPERTIES(${_NAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
      IF(WIN32 AND NOT CYGWIN)
        SET_TARGET_PROPERTIES(${_NAME} PROPERTIES SUFFIX ".pyd")
      ENDIF(WIN32 AND NOT CYGWIN)
    ENDIF(PYTHON_MODULE_${_NAME}_BUILD_SHARED)

  ENDIF(PYTHON_ENABLE_MODULE_${_NAME})
ENDFUNCTION(PYTHON_ADD_MODULE)

FUNCTION(PYTHON_WRITE_MODULES_HEADER _filename)

  GET_PROPERTY(PY_STATIC_MODULES_LIST  GLOBAL  PROPERTY PY_STATIC_MODULES_LIST)

  GET_FILENAME_COMPONENT(_name "${_filename}" NAME)
  STRING(REPLACE "." "_" _name "${_name}")
  STRING(TOUPPER ${_name} _nameUpper)
  SET(_filename ${CMAKE_CURRENT_BINARY_DIR}/${_filename})

  SET(_filenameTmp "${_filename}.in")
  FILE(WRITE ${_filenameTmp} "/*Created by cmake, do not edit, changes will be lost*/\n")
  FILE(APPEND ${_filenameTmp}
"#ifndef ${_nameUpper}
#define ${_nameUpper}

#include <Python.h>

#ifdef __cplusplus
extern \"C\" {
#endif /* __cplusplus */

")

  FOREACH(_currentModule ${PY_STATIC_MODULES_LIST})
    FILE(APPEND ${_filenameTmp} "extern void init${PYTHON_MODULE_PREFIX}${_currentModule}(void);\n\n")
  ENDFOREACH(_currentModule ${PY_STATIC_MODULES_LIST})

  FILE(APPEND ${_filenameTmp}
"#ifdef __cplusplus
}
#endif /* __cplusplus */

")


  FOREACH(_currentModule ${PY_STATIC_MODULES_LIST})
    FILE(APPEND ${_filenameTmp} "int ${_name}_${_currentModule}(void) \n{\n  static char name[]=\"${PYTHON_MODULE_PREFIX}${_currentModule}\"; return PyImport_AppendInittab(name, init${PYTHON_MODULE_PREFIX}${_currentModule});\n}\n\n")
  ENDFOREACH(_currentModule ${PY_STATIC_MODULES_LIST})

  FILE(APPEND ${_filenameTmp} "void ${_name}_LoadAllPythonModules(void)\n{\n")
  FOREACH(_currentModule ${PY_STATIC_MODULES_LIST})
    FILE(APPEND ${_filenameTmp} "  ${_name}_${_currentModule}();\n")
  ENDFOREACH(_currentModule ${PY_STATIC_MODULES_LIST})
  FILE(APPEND ${_filenameTmp} "}\n\n")
  FILE(APPEND ${_filenameTmp} "#ifndef EXCLUDE_LOAD_ALL_FUNCTION\nvoid CMakeLoadAllPythonModules(void)\n{\n  ${_name}_LoadAllPythonModules();\n}\n#endif\n\n#endif\n")

# with CONFIGURE_FILE() cmake complains that you may not use a file created using FILE(WRITE) as input file for CONFIGURE_FILE()
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_filenameTmp}" "${_filename}" OUTPUT_QUIET ERROR_QUIET)

ENDFUNCTION(PYTHON_WRITE_MODULES_HEADER)
