if(ENABLE_cuEST)
    find_package(cuEST 0.1 CONFIG QUIET)  # edit in codedeps

    if(NOT cuEST_FOUND)
        find_package(cuEST MODULE QUIET)
    endif()

    if(cuEST_FOUND)
        message(STATUS "${Cyan}Found cuEST${ColourReset}: ${cuEST_LIBRARY} (include: ${cuEST_INCLUDE_DIR})")
        add_library(cuest_external INTERFACE)  # dummy
    else()
        if(${CMAKE_INSIST_FIND_PACKAGE_cuEST})
            message(FATAL_ERROR "Suitable cuEST could not be externally located as user insists")
        endif()

        message(STATUS "Suitable cuEST could not be located, set -DcuEST_ROOT=/path/to/cuest or ENABLE_cuEST=OFF")
        add_library(cuest_external INTERFACE)  # dummy
    endif()
else()
    message(STATUS "Optional addon ${Blue}unsought: cuEST${ColourReset}")
    add_library(cuest_external INTERFACE)  # dummy
endif()
