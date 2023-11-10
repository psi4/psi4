find_package(BrianQC 1.1 CONFIG QUIET)

if(${BrianQC_FOUND})
    get_property(_loc TARGET BrianQC::static_wrapper PROPERTY LOCATION)
    message(STATUS "${Cyan}Found BrianQC${ColourReset}: ${_loc} (found version ${BrianQC_VERSION})")
    add_library(brianqc_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_BrianQC})
        message(VERBOSE "Suitable BrianQC located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_BrianQC})
        message(FATAL_ERROR "Suitable BrianQC could not be externally located as user insists")
    endif()

    message(STATUS "Suitable BrianQC could not be located, please contact info@brianqc.com to obtain BrianQC")

    add_library(brianqc_external INTERFACE)  # ensuring that the dummy target exists anyway
endif()
