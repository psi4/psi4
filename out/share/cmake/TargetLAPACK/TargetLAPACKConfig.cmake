
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was TargetLAPACKConfig.cmake.in                            ########

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

set(PN TargetLAPACK)
check_required_components(${PN})

# make detectable the various cmake modules exported alongside
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

set(isMKL " MKL")

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET tgt::lapack)
    include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets.cmake")

    include(CMakeFindDependencyMacro)
    if(NOT TARGET tgt::MathOpenMP)
        find_dependency(MathOpenMP)
    endif()

    get_property(_ill TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property(_illl TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property(_illb TARGET tgt::blas PROPERTY INTERFACE_LINK_LIBRARIES)

    #include(CMakePrintHelpers)
    #cmake_print_properties(TARGETS OpenMP::OpenMP_C OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran OpenMP::OpenMP tgt::MathOpenMP tgt::lapack
    #                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)
endif()
