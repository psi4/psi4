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

include_directories(include)
include_directories(${PSI4_ROOT})
add_subdirectory(src)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
add_executable(psi4_exe psi4.cc)

get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARIES} DIRECTORY)

set(psi4_install_lib ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(PSI4_RPATH ${psi4_install_lib}
        ${DESTDIR}/${psi4_install_lib}
        ${PYTHON_LIB_DIR}
        )
set_target_properties(psi4_exe PROPERTIES OUTPUT_NAME psi4
        INSTALL_RPATH "${PSI4_RPATH}"
        INSTALL_RPATH_USE_LINK_PATH TRUE
        BUILD_WITH_INSTALL_RPATH TRUE
        )
get_property(psi4_libraries GLOBAL PROPERTY LIBLIST)
message("LIBS: PYTHON   ${PYTHON_LIBRARIES}")
message("LIBS: PSI4     ${psi4_libraries}")
message("LIBS: CHEMPS2  ${CHEMPS2_LIBRARIES}")
message("LIBS: GDMA")
message("LIBS: DERIV    ${LIBDERIV_LIBRARIES}")
message("LIBS: THREAD   ${CMAKE_THREAD_LIBS_INIT}")
message("LIBS: LAPACK   ${LAPACK_LIBRARIES}")
message("LIBS: BLAS     ${BLAS_LIBRARIES}")
message("LIBS: Fortran  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

target_link_libraries(psi4_exe PRIVATE ${PYTHON_LIBRARIES})
target_link_libraries(psi4_exe PRIVATE psi4_main)
target_link_libraries(psi4_exe PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(psi4_exe PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(psi4_exe PRIVATE ${BLAS_LIBRARIES})

#get_property(psi4_modules GLOBAL PROPERTY PSI4_MODULES)
#why twicetarget_link_libraries(psi4_exe PRIVATE ${psi4_lib_list} psi4_main)
#add_dependencies(psi4_exe psi4_main)
#set(psi4_lib_list ${PYTHON_LIBRARIES}
#                  psi4_main
#                  ${psi4_modules}
#                  ${LIBDERIV_LIBRARIES}
#                  ${CMAKE_THREAD_LIBS_INIT}
#                  ${LAPACK_LIBRARIES}
#                  ${BLAS_LIBRARIES}
#)
#target_link_libraries(psi4_exe PRIVATE ${psi4_lib_list} psi4_main)

configure_file("psi4-config.in" "psi4-config.tmp" @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4-config
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

configure_file("psi4Config.cmake.in" "psi4Config.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/psi4Config.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)
install(DIRECTORY share DESTINATION ${CMAKE_INSTALL_PREFIX}
        MESSAGE_NEVER
        PATTERN fsapt EXCLUDE
        )
install(TARGETS psi4_exe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
