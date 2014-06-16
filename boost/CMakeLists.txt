add_custom_command(
    COMMAND tar xzf ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz 
    OUTPUT ${BOOSTVER}
    COMMENT Unpacking Boost
    DEPENDS ${PROJECT_SOURCE_DIR}/boost/${BOOSTVER}.tar.gz
)
set(MPI_SENT "")
if(MPI_FOUND)
  set(BOOST_LIBS_2_MAKE
"filesystem,mpi,python,regex,serialization,system,thread")
  set(MPI_SENT "using mpi \\;")
else()
  set(BOOST_LIBS_2_MAKE "filesystem,python,regex,serialization,system,thread")
endif()

#To get boost to compile MPI we need to append "using mpi ;" to the end of the
#project-config.jam file

add_custom_command(
    COMMAND ./bootstrap.sh --with-libraries=${BOOST_LIBS_2_MAKE}
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
