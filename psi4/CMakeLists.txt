cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
# no Fortran in psi4-core proper, but language needs to be declared
#   for CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to be populated so
#   Fortran add-ons can be linked
project(psi4-core LANGUAGES C CXX Fortran VERSION 1.0)

#set(CMAKE_CXX_STANDARD 11) enable if they fix this for Intel
list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)

include(Psi4Macros)
include(TestRestrict)
test_restrict(restrict)

### >> Locate Python Interpreter and Library << ###
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

### >> Locate Boost Components << ###
set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
        thread)
find_package(Boost 1.57.0 QUIET REQUIRED COMPONENTS ${BOOST_COMPONENTS})

set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})
find_package(Threads)
find_package(DL)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    find_package(RT REQUIRED)
    list(APPEND Boost_LIBRARIES ${LIBRT_LIBRARIES})
endif ()

#find_package(libderiv REQUIRED COMPONENTS ${LIBINT_OPT_AM})
find_package(libint CONFIG REQUIRED COMPONENTS ${MAX_AM_ERI}
        HINTS ${CMAKE_PREFIX_PATH}/libint)

set_property(GLOBAL PROPERTY PSI4_MODULES "")
set_property(GLOBAL PROPERTY LIBLIST "")

#optional_plugin(Ambit)
#optional_plugin(PCMSolver)
#optional_plugin(CHEMPS2)
#optional_plugin(liberd)
#optional_plugin(dkh)
#optional_plugin(gdma)
#include(get_gdma)


if (${ENABLE_CHEMPS2})
    find_package(CHEMPS2 REQUIRED)
    #   set_property(GLOBAL APPEND PROPERTY PSI4_MODULES ${${PLUGIN_NAME}_LIBRARIES})
    add_definitions(-DENABLE_CHEMPS2)
endif ()

if (${ENABLE_gdma})
    find_package(gdma CONFIG REQUIRED
            HINTS ${CMAKE_PREFIX_PATH}/gdma)
    message("gdma found w/i ${CMAKE_PREFIX_PATH}")
else ()
    message("gdma not enabled ${CMAKE_PREFIX_PATH}")
endif ()
include(echo_target)
if (TARGET gdma::gdma)
    print_target_vars(gdma)
    print_target(gdma::gdma)
endif ()
find_package(libefp CONFIG REQUIRED
        HINTS ${CMAKE_PREFIX_PATH}/libefp)
message("libefp found w/i ${CMAKE_PREFIX_PATH}")
print_target_vars(libefp)
print_target(libefp::efp)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(include)
include_directories(src)
include_directories(${PSI4_ROOT})

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(include)
add_subdirectory(src)

configure_file("psi4-config.in" "psi4-config.tmp" @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4-config
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file("psi4Config.cmake.in" "psi4Config.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/psi4Config.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY share DESTINATION ${CMAKE_INSTALL_PREFIX}
        MESSAGE_NEVER
        PATTERN fsapt EXCLUDE
        )
