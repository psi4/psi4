project(TargetLAPACK)

# <<<  "Build"  >>>

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/math)

set(MKL_COMPILER_BINDINGS "${CMAKE_CXX_COMPILER_ID}")
set(BLAS_LANG "CXX")
set(LAPACK_LANG "CXX")
if(NOT LAPACK_LIBRARIES)
    include(ConfigMath)
endif()

# <<  Build targets  >>
add_library(blas INTERFACE)
add_library(lapk INTERFACE)
add_library(lapack INTERFACE)
set_property(TARGET blas PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
set_property(TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
set_property(TARGET blas PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS})
set_property(TARGET lapk PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${LAPACK_INCLUDE_DIRS})
set_property(TARGET lapack PROPERTY INTERFACE_LINK_LIBRARIES lapk blas tgt::MathOpenMP)

# <<  Examine library list for MKL-ness  >>
foreach(_l IN LISTS LAPACK_LIBRARIES)
    get_filename_component(_lname ${_l} NAME)
    if(${_lname} MATCHES "mkl")
        set(isMKL " MKL")
        target_compile_definitions(lapack INTERFACE USING_LAPACK_MKL)
        break()
    endif()
endforeach()
if(isMKL STREQUAL " MKL")
    set(_VENDOR "MKL")
else()
    set(_VENDOR "All")
endif()
set_target_properties(
  lapack
  PROPERTIES
    VENDOR ${_VENDOR}
    INT_INTERFACE lp64  # Psi4 only works with lp64 afaik
  )
set(export_properties "VENDOR;INT_INTERFACE")
set_property(TARGET lapack APPEND PROPERTY EXPORT_PROPERTIES "${export_properties}")

# <<  Detect OpenMP and modify for BLAS/LAPACK  >>
if(NOT TARGET tgt::MathOpenMP)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
    set(TargetOpenMP_FIND_COMPONENTS "CXX")
    find_package(MathOpenMP)
endif()

#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS OpenMP::OpenMP_CXX OpenMP::OpenMP tgt::MathOpenMP lapack
#                       PROPERTIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRS INTERFACE_LINK_LIBRARIES)

get_property(_ill TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES)
list(GET _ill 0 _ill0)
message(STATUS "${Cyan}Found LAPACK${isMKL}${ColourReset}: ${_ill0};...")

if(NOT (${isMKL} MATCHES "MKL"))
    message(WARNING "${Yellow}MKL is the only BLAS/LAPACK distribution thoroughly tested with Psi4. Several others are mostly fine, with problems cropping up in particular modules or when running in threaded mode. When using a non-MKL distribution, please run the test suite in full and also run particular tests in threaded mode, if you plan to run Psi4 threaded. MKL is freely available through conda; see instructions to compile against conda MKL through the 'source' option on the download website page. OpenBLAS also seems to work, but be sure you're using the OpenMP, not the pthreads, variant.${ColourReset}")
endif()

unset(BLAS_LIBRARIES)
unset(LAPACK_LIBRARIES)
unset(BLAS_INCLUDE_DIRS)
unset(LAPACK_INCLUDE_DIRS)

set(PN ${PROJECT_NAME})
install(TARGETS blas lapk lapack
        EXPORT "${PN}Targets")

# <<<  Export Config  >>>

include(CMakePackageConfigHelpers)

set(CMAKECONFIG_INSTALL_DIR "${TargetLAPACK_INSTALL_CMAKEDIR}")
configure_package_config_file(${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_SOURCE_DIR}/FindTargetOpenMP.cmake
              ${CMAKE_CURRENT_SOURCE_DIR}/FindMathOpenMP.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "tgt::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

add_custom_target(lapack_external
    ${CMAKE_COMMAND} "-DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}" -P "${PROJECT_BINARY_DIR}/cmake_install.cmake"
    COMMENT "Installing Psi4-detected BLAS/LAPACK")

set(${PN}_DIR ${STAGED_INSTALL_PREFIX}/${TargetLAPACK_INSTALL_CMAKEDIR} CACHE PATH "path to externally detected ${PN}Config.cmake" FORCE)
