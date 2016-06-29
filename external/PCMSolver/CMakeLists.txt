include(ExternalProject)

find_package(PCMSolver)
if(NOT PCMSolver_FOUND AND ENABLE_PCMSOLVER)
  message(STATUS "PCMSolver not found. The pre-packaged version will be built.")

  #find_package(ZLIB QUIET)
  #if (NOT ZLIB_FOUND)
  #  message(FATAL_ERROR "No Zlib, no PCMSolver. Build against existing with -DPCMSOLVER_ROOT=/path/to/pcmsolver or skip with -DENABLE_PCMSOLVER=OFF")
  #endif()

  #set(CUSTOM_PCMSolver_LOCATION ${PROJECT_BINARY_DIR}/interfaces/pcmsolver)

  # Make sure PCMSolver gets exact same Boost as Psi4
  if(BUILD_CUSTOM_BOOST)
    ExternalProject_Add(pcmsolver_external
       DEPENDS boost_external
       GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
       GIT_TAG v1.1.1
       CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pcmsolver
                  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                  -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
                  -DENABLE_CXX11_SUPPORT=ON
                 #-DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
                 #-DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
                  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
                  -DENABLE_TESTS=OFF
                  -DENABLE_LOGGER=OFF
                  -DENABLE_TIMER=OFF
                  -DBUILD_STANDALONE=OFF
                  -DENABLE_FORTRAN_API=OFF
                  -DSTATIC_LIBRARY_ONLY=ON
                  -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
                 #-DZLIB_ROOT=${ZLIB_ROOT}
                  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABL
       INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                      DESTDIR=${CMAKE_BINARY_DIR}/stage
    )
else()
   add_library(pcmsolver_external INTERFACE)
endif()

# [1] It's nice to have a full PCMSolver::PCMSolver target that has embedded
# the library, the linking library paths with dependencies, and the include
# paths with dependencies, just like FindPCMSolver supplies, especially
# since src/bin/libpsipcm needs the headers for compilation and src/bin/psi4
# needs all the linking libraries. Problem is that conventional target
# derived from ExternalProject_Add is a highly sought but not quite
# certified cmake pattern. Hence INTERFACE_INCLUDE_DIRECTORIES complains
# that the directories don't exist at configure time. Hence the hack to
# create an empty directory.

