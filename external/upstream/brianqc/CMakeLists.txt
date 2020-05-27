find_package(BrianQC 1.1 CONFIG QUIET)

if(${BrianQC_FOUND})
    get_property(_loc TARGET BrianQC::static_wrapper PROPERTY LOCATION)
    if(${BrianQC_DEV})
        message(STATUS "${Cyan}Found BrianQC${ColourReset}: ${_loc} (found development version ${BrianQC_VERSION})")
    else(${BrianQC_DEV})
        message(STATUS "${Cyan}Found BrianQC${ColourReset}: ${_loc} (found installed version ${BrianQC_VERSION})")
    endif(${BrianQC_DEV})
    
    add_library(brianqc_external INTERFACE)  # dummy
else()
    if(${CMAKE_INSIST_FIND_PACKAGE_BrianQC})
        message(FATAL_ERROR "Suitable BrianQC could not be externally located as user insists")
    endif()
    
    message(STATUS "Suitable BrianQC could not be located, please contact info@brianqc.com to obtain BrianQC")
endif()
