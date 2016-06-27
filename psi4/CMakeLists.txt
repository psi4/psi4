cmake_minimum_required(VERSION 3.0)
project(psi4 LANGUAGES C CXX)
include(../cmake/Psi4Macros.cmake)
#Find boost doesn't appear to follow the usual search rules and will miss our
#compiled boost
set(BOOST_ROOT ${CMAKE_PREFIX_PATH}/boost)
set(BOOST_COMPONENTS filesystem python regex serialization system timer chrono
                     thread
)
find_package(Boost 1.57.0 QUIET REQUIRED COMPONENTS ${BOOST_COMPONENTS})
#For some stupid reason boost doesn't add its dependencies to it's list of 
#requirements like a good CMake program...
set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})
#et(Boost_LIBRARIES ${Boost_LIBRARIES} ${PYTHON_LIBRARY})
include_directories(${PSI4_ROOT})
add_subdirectory(src)
