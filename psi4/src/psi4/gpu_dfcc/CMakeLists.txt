cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string
project(gpu_dfcc 
        VERSION 0.1
        LANGUAGES CXX CUDA)
set(gpudfcc_AUTHORS      "A. Eugene DePrince III")
set(gpudfcc_DESCRIPTION  "GPU DF-CC")
#set(gpudfcc_URL          "github.com/edeprince3/gpu_dfcc")
set(gpudfcc_LICENSE      "GPL 2+")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# <<<  Options  >>>

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler --compiler-options '-fPIC' -O2 ${CUDA_LIBRARIES} -lcublas -lcudart -lcusolver")
#set(CXX_FLAGS ${CUDA_LIBRARIES} "-lcublas -lcudar -lcusolver")
include(psi4OptionsTools)
#set(CMAKE_CXX_STANDARD_LIBRARIES ${CUDA_LIBRARIES} "-lcublas -lcudart -lcusolver")
option_with_print(ENABLE_OPENMP "Enable threadsafe linking to OpenMP parallelized programs." ON)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
if (APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES GNU))
    option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" OFF)
else ()
    option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                      "-xHost" "-march=native")
endif()
#find_package(CUDA)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})
#
# PYMOD_INSTALL_LIBDIR If set, installs python modules to
#                      ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/v2rdm_casscf rather than
#                      the default read off parent Psi4

#   install alongside psi4 module by default, but overrideable
#get_filename_component(psi4_CMAKE_INSTALL_PREFIX ${psi4_INCLUDE_DIR} DIRECTORY)
#if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    set(CMAKE_INSTALL_PREFIX ${psi4_CMAKE_INSTALL_PREFIX} CACHE PATH "Install path" FORCE)
#endif()
#message(STATUS "Plugin install prefix: ${CMAKE_INSTALL_PREFIX}")
#
#if(PYMOD_INSTALL_LIBDIR)
#    set(PYMOD_INSTALL_FULLDIR "${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/${PN}")
#else()
#    file(RELATIVE_PATH _tmp ${psi4_CMAKE_INSTALL_PREFIX} ${psi4_LIBRARY})
    #   e.g., _tmp = lib/psi4/core.so
#    get_filename_component(_tmp2 ${_tmp} DIRECTORY)
#    get_filename_component(_tmp3 ${_tmp2} DIRECTORY)
#    set(PYMOD_INSTALL_FULLDIR "${_tmp3}/${PN}")
#endif()
#message(STATUS "Plugin module install: ${PYMOD_INSTALL_FULLDIR}")


#include(custom_cxxstandard)
#include(autocmake_omp)
#include(custom_static_library)

# <<<  Build  >>>

#if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays")
#endif()
#set(CMAKE_CUDA_FLAGS "--std=c++11")
# old way of listing files
#file(GLOB sources "*.cc" "*.F90" "*.cu")
#cuda_add_library(gpu_dfcc 
#        ${SOURCES})
set (CMAKE_CXX_FLAGS "--std=c++11")
set (CUDA_PROPAGATE_HOST_FLAGS ON)
set(sources ccsd.cu gpu_dfcc.cu gpuhelper_driver.cc triples.cu)
set(libname gpu_dfcc)
message("${libname}")
    set(current_sources ${${sources}};)
    add_library(${libname} STATIC ccsd.cu gpu_dfcc.cu gpuhelper_driver.cc triples.cu)
    set_target_properties(${libname} PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})

    # library modules get their headers installed
    if(${binlib} MATCHES lib)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4
                FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.i")
    endif()

    # binary modules explicitly compiled into psi4.so
    if(${binlib} MATCHES bin)
        set_property(GLOBAL APPEND PROPERTY BINLIST ${libname})
    endif()

    set(depend_name "${ARGN}")
    foreach(name_i IN LISTS depend_name)
        target_link_libraries(${libname} ${name_i})
    endforeach()
    target_link_libraries(${libname} PRIVATE pybind11::module)
    target_link_libraries(${libname} PUBLIC tgt::lapack)
    target_link_libraries(${libname} PRIVATE core)
#add_psi4_plugin(gpu_dfcc
#)

#target_link_libraries(gpu_dfcc PRIVATE ${LIBC_INTERJECT})

# <<<  Install  >>>

#install(TARGETS gpu_dfcc
#        EXPORT "${PN}Targets"
#        LIBRARY DESTINATION ${PYMOD_INSTALL_FULLDIR})

#install(FILES __init__.py pymodule.py LICENSE README.md
#        DESTINATION ${PYMOD_INSTALL_FULLDIR})

#install(DIRECTORY tests/
#        DESTINATION ${PYMOD_INSTALL_FULLDIR}/tests
#        FILES_MATCHING PATTERN "gpu_dfcc*/input.dat"
#        PATTERN "old-tests" EXCLUDE)

# <<< Export Config >>>




# old stuff

#target_compile_definitions(gpu_dfcc PRIVATE _OPENMP)

