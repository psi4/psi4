# Just change the Boost version number here
set(BOOSTVER 1.57.0)
set(BOOSTVERMIN 1.55.0)
if (NOT DEFINED BUILD_CUSTOM_BOOST)
    set(BUILD_CUSTOM_BOOST FALSE)
endif()
# List all components needed (unit_test_framework) here.
# unit_test_framework will be added afterwards, if needed.
# RMR mpi has to be added here or else it's not added to the linking stage
list(APPEND needed_components filesystem python regex serialization system timer chrono thread)
if(ENABLE_MPI)
   list(APPEND needed_components mpi)
endif()
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
if(ENABLE_UNIT_TESTS)
   list(APPEND needed_components unit_test_framework)
endif()

if (NOT BUILD_CUSTOM_BOOST)
    find_package(Boost ${BOOSTVERMIN} COMPONENTS "${needed_components}")
else()
    set(Boost_FOUND FALSE)
    set(BUILD_CUSTOM_BOOST TRUE)
endif()
if(NOT Boost_FOUND)
   # Set also variables usually set by find_package
   message(STATUS "Boost ${BOOSTVERMIN} not found. The pre-packaged version will be built.")
   set(BUILD_CUSTOM_BOOST TRUE)
   set(CUSTOM_BOOST_LOCATION ${PROJECT_BINARY_DIR}/boost)
   string(REGEX REPLACE "\\." "0" Boost_VERSION ${BOOSTVER})
   math(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
   math(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
   math(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
   set(Boost_LIB_VERSION ${Boost_MAJOR_VERSION}_${Boost_MINOR_VERSION})
   add_subdirectory(boost)
   set(Boost_FOUND TRUE)
   set(Boost_LIBRARIES "")
   # Read documentation in FindBoost.cmake for the difference between the singular and plural forms
   set(Boost_INCLUDE_DIR  ${CUSTOM_BOOST_LOCATION}/include)
   set(Boost_INCLUDE_DIRS ${CUSTOM_BOOST_LOCATION}/include)
   set(Boost_LIBRARY_DIR  ${CUSTOM_BOOST_LOCATION}/lib)
   set(Boost_LIBRARY_DIRS ${CUSTOM_BOOST_LOCATION}/lib)
   # This is the one that was in use in the PSI4 cmake files and is included
   # to maintain the scripts in some working order, but the former should
   # be preferred
   set(BOOSTLIBDIR ${Boost_LIBRARY_DIR})
   # We will link statically, so just set the Boost_<C>_LIBRARY for the static library
   foreach(_component ${needed_components})
   string(TOUPPER ${_component} _COMP)
   set(Boost_${_COMP}_FOUND TRUE)
   set(Boost_${_COMP}_LIBRARY libboost_${_component}-${Boost_LIB_VERSION}.a)
   list(APPEND Boost_LIBRARIES ${Boost_${_COMP}_LIBRARY})
   endforeach()
   if(CMAKE_SYSTEM_NAME MATCHES "Linux")
      list(APPEND Boost_LIBRARIES rt)
   endif()
endif()
