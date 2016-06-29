cmake_minimum_required(VERSION 3.0)
project(psi4 LANGUAGES C CXX VERSION 1.0)
list(APPEND CMAKE_MODULE_PATH ${PSI4_ROOT}/cmake)
include(Psi4Macros)

set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
                     thread
)
find_package(Boost 1.57.0 QUIET REQUIRED COMPONENTS ${BOOST_COMPONENTS})
#For some stupid reason boost doesn't add its dependencies to it's list of 
#requirements like a good CMake program...
set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})

optional_plugin(Ambit)
optional_plugin(PCMSolver)
optional_plugin(CHEMPS2)

include_directories(${PSI4_ROOT})
add_subdirectory(src)

#add_executable(psi4 psi4.cc)
#target_link_libraries(psi4 ${LIBLIST} ${PSI4_MODULES})