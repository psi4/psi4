# FindTargetOpenMP.cmake
# ----------------------
#
# OpenMP cmake module to wrap FindOpenMP.cmake in a target.
#
# This module sets the following variables in your project: ::
#
#   TargetOpenMP_FOUND - true if OpenMP found on the system
#   TargetOpenMP_MESSAGE - status message with OpenMP library path list
#
# * this is a hack until we update min cmake (to 3.10?)
# * it is specific to Psi4 in that:
#   * we only care about C++ OpenMP
#   * minimum Intel in 2017, so no need to add logic for lower
#   * we never use OPENMP_FOUND
# * this is run by the TargetLAPACK module and tgt::OpenMP is supplemented by
#   math libraries detected there

set(PN TargetOpenMP)

# 1st precedence - libraries passed in through -DOpenMP_LIBRARIES
if (OpenMP_LIBRARIES AND OpenMP_FLAGS)
    if (NOT ${PN}_FIND_QUIETLY)
        message (STATUS "OpenMP detection suppressed.")
    endif()

    add_library (OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_FLAGS})
    set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_LIBRARIES})
else()
    # 2nd precedence - target from FindOpenMP.cmake
    find_package (OpenMP QUIET MODULE COMPONENTS CXX)
    if (TARGET OpenMP::OpenMP_CXX)
    else()
        # 3rd precedence - construct a target
        add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)

        if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:-fopenmp>")
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES "gomp;pthread")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:-qopenmp>")
            set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES "iomp5;pthread")
        endif()
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "OpenMP target constructed.")
        endif()
    endif()
endif()

get_property (_ill TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES)
set (${PN}_MESSAGE "Found TargetOpenMP: ${_ill}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (${PN} DEFAULT_MSG ${PN}_MESSAGE)
