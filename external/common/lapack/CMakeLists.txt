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
add_library(omp INTERFACE)
add_library(blas INTERFACE)
add_library(lapk INTERFACE)
add_library(lapack INTERFACE)
set_property(TARGET omp PROPERTY INTERFACE_LINK_LIBRARIES ${OMP_LIBRARIES})
set_property(TARGET blas PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
set_property(TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
set_property(TARGET omp PROPERTY INTERFACE_COMPILE_OPTIONS ${OMP_FLAGS})
set_property(TARGET blas PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS})
set_property(TARGET lapk PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${LAPACK_INCLUDE_DIRS})
set_property(TARGET lapack PROPERTY INTERFACE_LINK_LIBRARIES lapk blas omp)
get_property(_ill TARGET lapk PROPERTY INTERFACE_LINK_LIBRARIES)

if(ENABLE_OPENMP)
    if(TARGET OpenMP::OpenMP_CXX)
    else()
        list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
        find_package(TargetOpenMP COMPONENTS CXX)
        target_link_libraries(omp INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

foreach(_l IN LISTS LAPACK_LIBRARIES)
    if(${_l} MATCHES "mkl")
        set(_isMKL " MKL")
        target_compile_definitions(lapack INTERFACE -DUSING_LAPACK_MKL)
        break()
    endif()
endforeach()

# <<  Examine library list for MKL-ness  >>
list(GET _ill 0 _ill0)
message(STATUS "${Cyan}Found LAPACK${_isMKL}${ColourReset}: ${_ill0};...")

unset(BLAS_LIBRARIES)
unset(LAPACK_LIBRARIES)
unset(BLAS_INCLUDE_DIRS)
unset(LAPACK_INCLUDE_DIRS)

set(PN ${PROJECT_NAME})
install(TARGETS omp blas lapk lapack
        EXPORT "${PN}Targets")

# <<<  Export Config  >>>

include(CMakePackageConfigHelpers)

set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_SOURCE_DIR}/FindTargetOpenMP.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "tgt::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

add_custom_target(lapack_external
    ${CMAKE_COMMAND} "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}" -P "${PROJECT_BINARY_DIR}/cmake_install.cmake"
    COMMENT "Installing Psi4-detected BLAS/LAPACK")

set(${PN}_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/${PN} CACHE PATH "path to externally detected ${PN}Config.cmake" FORCE)
