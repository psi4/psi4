project(psi4)
cmake_minimum_required(VERSION 2.8)


#
# Python Detection 
#
find_package(PythonInterp REQUIRED)
set(PYCONFIG ${PYTHON_EXECUTABLE}-config)
execute_process(COMMAND ${PYCONFIG} --libs     OUTPUT_VARIABLE PYLIB OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PYCONFIG} --include  OUTPUT_VARIABLE PYINC OUTPUT_STRIP_TRAILING_WHITESPACE)
if(PYLIB STREQUAL "" OR PYINC STREQUAL "")
   message(FATAL_ERROR "Python development libraries were not detected. Make sure that the development version "
                       "is installed and that python-config exists in the same directory as the python executable.")
endif()

#
# Boost Detection 
#
# N.B. Update boost/CMakeLists.txt if this list of components changes
set(BOOST_COMPONENTS filesystem python regex thread system serialization)
find_package(Boost 1.47.0 COMPONENTS ${BOOST_COMPONENTS})
if(Boost_FOUND)
    # Use the Boost installation detected on this system
    set(BOOSTLIB ${Boost_LIBRARIES})
    set(BOOSTINC ${Boost_INCLUDE_DIRS})
else()
    # Add the Boost version bundled with Psi4 as a dependency
    message(WARNING "Boost not found. The pre-packaged version will be built.")
    # The version bundled with Psi4
    set(BOOSTVER boost_1_53_0)
    # Link and include information
    set(BOOSTLIBDIR ${PROJECT_BINARY_DIR}/boost/${BOOSTVER}/stage/lib)
    set(BOOSTINC ${PROJECT_BINARY_DIR}/boost)
    foreach(component ${BOOST_COMPONENTS})
        set(BOOSTLIB ${BOOSTLIB} ${BOOSTLIBDIR}/${component})
    endforeach()
    add_subdirectory(boost)
endif()

#
# Libint configuration
#
set(LIBINT_TMPDIR /tmp/libint)
set(LIBINT_OPT_AM 5)
math(EXPR LIBINT_NEW_AM ${LIBINT_OPT_AM}*2)
configure_file(src/lib/libint/libint_config.h.in src/lib/libint/libint_config.h)

# Psi4-specific includes and libraries
set(PSIINC ${PROJECT_BINARY_DIR}/include 
           ${PROJECT_SOURCE_DIR}/include
           ${PROJECT_SOURCE_DIR}/src/lib
           ${PROJECT_BINARY_DIR}/src/lib)

set(CMAKE_INCLUDE_CURRENT_DIR true)
include_directories(${BOOSTINC} ${PSIINC} ${PYINC})

add_subdirectory(src)


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
