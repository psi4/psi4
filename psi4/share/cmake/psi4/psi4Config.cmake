# psi4Config.cmake
# ----------------
#
# Psi4 cmake module.
# This module sets the following variables in your project::
#
#   psi4_FOUND - true if psi4 and all required components found on the system
##   psi4_VERSION - psi4 version in format Major.Minor.Release
##   psi4_INCLUDE_DIRS - Directory where psi4 header is located.
#   psi4_INCLUDE_DIR - Directory where psi4/psi4-dec.h is located.
##   psi4_DEFINITIONS - Definitions necessary to use psi4, namely USING_psi4.
##   psi4_LIBRARIES - psi4 library to link against.
#   psi4_LIBRARY - core.cpython-312-x86_64-linux-gnu.so library to link against.
#   psi4_EXECUTABLE - psi4 executable script.
##   psi4_FRAGLIB_DIRS - Directories (list) where EFP fragments are located
#
#
# Available components::
#
#   ambit - tensor library active
#   chemps2 - dmrg library active
#   dkh - relativistic library active
#   ecpint - ecp integrals library active
#   gdma - multipole library active
#   pcmsolver - solvation library active
#   simint - integrals library active
#   einsums - tensor contraction library active
#   gauxc - exact exchange library active 
#
#
# Exported targets::
#
# If psi4 is found, this module defines the following :prop_tgt:`IMPORTED`
# target. ::
#
#   psi4::core - the main psi4 library with header & pybind11 & Python headers attached.
#
#
# Suggested usage::
#
#   find_package(psi4)
##   find_package(psi4 1.3.0 EXACT CONFIG REQUIRED COMPONENTS shared)
#
#
# The following variables can be set to guide the search for this package::
#
#   psi4_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   PATH - environment variable, set to bin directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_psi4 - CMake variable, disables
#     find_package(psi4) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was psi4Config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(PN psi4)
set (_valid_components
    ambit
    chemps2
    dkh
    ecpint
    gdma
    pcmsolver
    simint
    einsums
    gauxc
)

# location of script psi4
set(${PN}_EXECUTABLE "${PACKAGE_PREFIX_DIR}/bin/psi4")

# directory containing psi4/psi4-dec.h
set(${PN}_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")

# location of psi4/core.so
set(${PN}_LIBRARY "${PACKAGE_PREFIX_DIR}/lib//psi4/core.cpython-312-x86_64-linux-gnu.so")
set(${PN}_SYS_PATH "${PACKAGE_PREFIX_DIR}/lib/")

# presence of external addons
macro(psi4_component name on)
    if(${on})
        set(${PN}_${name}_FOUND 1)
        list(APPEND ${PN}_FOUND_COMPONENTS ${name})
    endif()
endmacro()

psi4_component(ambit OFF)
psi4_component(chemps2 OFF)
psi4_component(dkh OFF)
psi4_component(ecpint ON)
psi4_component(gdma OFF)
psi4_component(pcmsolver OFF)
psi4_component(simint OFF)
psi4_component(einsums OFF)
psi4_component(gauxc OFF)

foreach(_comp ${${PN}_FIND_COMPONENTS})
    if(${${PN}_${_comp}_FOUND})
    else()
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN} missing requested component: ${_comp}")
        endif()
    endif()
endforeach()

check_required_components(${PN})

if(NOT CMAKE_REQUIRED_QUIET)
    message(STATUS "Psi4 script:   ${${PN}_EXECUTABLE}")
    message(STATUS "Psi4 headers:  ${${PN}_INCLUDE_DIR}")
    message(STATUS "Psi4 library:  ${${PN}_LIBRARY}")
    message(STATUS "Psi4 sys.path: ${${PN}_SYS_PATH}")
    message(STATUS "Psi4 components: ${${PN}_FOUND_COMPONENTS}")
    message(STATUS "Python executable: /storage/home/hcoda1/9/jpederson6/.conda/envs/qmmm/bin/python3.12")
endif()

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${PN}::core)
    include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets.cmake")
    set(CMAKE_CXX_STANDARD 20)
    if (${CMAKE_MINIMUM_REQUIRED_VERSION} VERSION_LESS 3.15)
        message(WARNING "Calling program/plugin using old Python detection because cmake_minimum_required (${CMAKE_MINIMUM_REQUIRED_VERSION}) is <3.15. Contact the program/plugin developer to update.")
    else ()
        find_package(Python 3.6 COMPONENTS Interpreter Development NumPy REQUIRED)
    endif()
    find_package(pybind11 CONFIG REQUIRED)
    if(NOT DEFINED ENABLE_OPENMP)
      set(ENABLE_OPENMP "ON")
    endif()
    find_package(TargetLAPACK CONFIG REQUIRED)
    find_package(Libint2 CONFIG REQUIRED)
    # No more ${PN} since defined in above find_package
endif()

function(add_psi4_plugin TARGET SOURCES)
    # remove ${TARGET} from ${ARGV} to use ${ARGV} as SOURCES
    list(REMOVE_AT ARGV 0)

    set(psi4_CXX_STANDARD 20)
    # checks compiler (and gcc, if necessary) C++14 compliant
    include(custom_cxxstandard)

    add_library(${TARGET} MODULE ${ARGV})
    target_link_libraries(${TARGET} PRIVATE psi4::core)
    target_link_libraries(${TARGET} PRIVATE tgt::MathOpenMP)
    target_include_directories(${TARGET} PRIVATE "${Python_INCLUDE_DIRS}")
    target_include_directories(${TARGET} PRIVATE "${PYTHON_INCLUDE_DIRS}")
    set_target_properties(${TARGET} PROPERTIES PREFIX "")
    if(APPLE)
      set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
endfunction()

# make detectable the various cmake modules exported by psi4
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
