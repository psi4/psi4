add_custom_command(
    COMMAND tar -xzf ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz 
    OUTPUT ${BOOSTVER}
    COMMENT Unpacking Boost
    DEPENDS ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz
)

#To get boost to compile MPI we need to append "using mpi ;" to the end of the
#user-config.jam file.  MPI_SENT will be the command we append
set(MPI_SENT "")
if(MPI_FOUND)
  set(BOOST_LIBS_2_MAKE
"filesystem,mpi,python,regex,serialization,system,thread")
  set(MPI_SENT "using mpi : ${CMAKE_CXX_COMPILER} \\;")
else()
  set(BOOST_LIBS_2_MAKE "filesystem,python,regex,serialization,system,thread")
endif()

#To get boost to be built with our compiler of choice we need these logic
#statements, then we write the options to user-config.jam
set(toolset "")
if(${CMAKE_CXX_COMPILER} MATCHES "icc*")
  if(UNIX OR APPLE)
    set(toolset "intel-linux")
  else()
    set(toolset "intel-windows")
  endif()
else()
   set(toolset "gcc")
endif()

add_custom_command(
    COMMAND echo "using ${toolset} : : ${CMAKE_CXX_COMPILER} \\;">user-config.jam
    COMMAND echo "${MPI_SENT}">>user-config.jam
    COMMAND ./bootstrap.sh --with-toolset=${toolset}
            --with-libraries=${BOOST_LIBS_2_MAKE}
           --with-python=${PYTHON}
    COMMAND echo "${MPI_SENT}">>project-config.jam
    COMMAND ./b2
    COMMAND touch ../boost.compiled
    WORKING_DIRECTORY ${BOOSTVER}
    DEPENDS ${BOOSTVER}
    OUTPUT boost.compiled
    OUTPUT ${BOOSTLIB}
    COMMENT Configuring and building Boost 
)

# This target is always executed, but checks for the existence of boost.compiled
add_custom_target(
    boostlibs ALL
    DEPENDS boost.compiled
)

#include(ExternalProject)
#ExternalProject_Add(Boost
#    URL ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz
#    CONFIGURE_COMMAND ./bootstrap.sh --with-libraries=${BOOST_LIBS_2_MAKE}
#    BUILD_COMMAND b2)
