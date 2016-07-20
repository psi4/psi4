cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(psi4 LANGUAGES C CXX VERSION 1.0)

#set(CMAKE_CXX_STANDARD 11) enable if they fix this for Intel
list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)
include(Psi4Macros)
include(TestRestrict)
test_restrict(restrict)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
                     thread
)
find_package(Boost 1.57.0 QUIET REQUIRED COMPONENTS ${BOOST_COMPONENTS})
#For some stupid reason boost doesn't add its dependencies to it's list of 
#requirements like a good CMake program...
set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})
find_package(Threads)
find_package(DL)
find_package(LAPACK REQUIRED)
find_package(RT REQUIRED)
list(APPEND Boost_LIBRARIES ${LIBRT_LIBRARIES})#again boost has failed...
find_package(libderiv REQUIRED COMPONENTS ${LIBINT_OPT_AM})
set_property(GLOBAL PROPERTY PSI4_MODULES "")
set_property(GLOBAL PROPERTY LIBLIST "")

optional_plugin(Ambit)
optional_plugin(PCMSolver)
optional_plugin(CHEMPS2)
optional_plugin(liberd)
optional_plugin(dkh)

include_directories(${PSI4_ROOT})
add_subdirectory(src)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(psi4_exe psi4.cc)
get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARIES} DIRECTORY)
set(PSI4_RPATH ${CMAKE_INSTALL_PREFIX}/lib
               ${PYTHON_LIB_DIR}
)
set_target_properties(psi4_exe PROPERTIES OUTPUT_NAME psi4
                          INSTALL_RPATH "${PSI4_RPATH}"
)
get_property(psi4_libraries GLOBAL PROPERTY LIBLIST)
get_property(psi4_modules GLOBAL PROPERTY PSI4_MODULES)
add_dependencies(psi4_exe psi4_main)
set(psi4_lib_list ${PYTHON_LIBRARIES}
                  psi4_main 
                  ${psi4_modules}
                  ${LIBDERIV_LIBRARIES}
                  ${CMAKE_THREAD_LIBS_INIT}
                  ${LAPACK_LIBRARIES}
)
target_link_libraries(psi4_exe PRIVATE ${psi4_lib_list})   


configure_file("psi4Config.cmake.in" "psi4Config.cmake" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/psi4Config.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX}/psi4)
install(DIRECTORY share DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS psi4_exe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
