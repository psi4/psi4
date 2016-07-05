
# User signal to try pre-built PCMSolver
if(PCMSOLVER_ROOT)
  find_package(PCMSolver)
endif()

# Build PCMSolver as external package if pre-built failed or not signaled
if(NOT PCMSolver_FOUND)
  message(STATUS "PCMSolver not found. The pre-packaged version will be built.")

  find_package(ZLIB QUIET)
  if (NOT ZLIB_FOUND)
    message(FATAL_ERROR "No Zlib, no PCMSolver. Build against existing with -DPCMSOLVER_ROOT=/path/to/pcmsolver or skip with -DENABLE_PCMSOLVER=OFF")
  endif()

  set(CUSTOM_PCMSolver_LOCATION ${PROJECT_BINARY_DIR}/interfaces/pcmsolver)
  include(ExternalProject)
  get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)
  # PCMSolver does not know profile
  if(CMAKE_BUILD_TYPE MATCHES "profile")
    set(PCM_BUILD_TYPE "release")
  else()
    set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  endif()
  list(APPEND PCMSolverCMakeArgs
    -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/interfaces
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DEXTRA_FCFLAGS=${PCM_EXTRA_Fortran_FLAGS}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DEXTRA_CFLAGS=${PCM_EXTRA_C_FLAGS}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DEXTRA_CXXFLAGS=${PCM_EXTRA_CXX_FLAGS}
    -DENABLE_CXX11_SUPPORT=${ENABLE_CXX11_SUPPORT}
    -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS}
    -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    -DENABLE_TESTS=OFF
    -DENABLE_LOGGER=OFF
    -DENABLE_TIMER=OFF
    -DBUILD_STANDALONE=OFF
    -DENABLE_FORTRAN_API=OFF
    -DSTATIC_LIBRARY_ONLY=ON
    -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
    -DZLIB_ROOT=${ZLIB_ROOT}
    -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    )
  # Make sure PCMSolver gets exact same Boost as Psi4
  if(BUILD_CUSTOM_BOOST)
    ExternalProject_Add(interface_pcmsolver
      DEPENDS custom_boost
      PREFIX ${CUSTOM_PCMSolver_LOCATION}
      GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
      GIT_TAG v1.1.3
      CMAKE_ARGS "${PCMSolverCMakeArgs}"
      INSTALL_DIR "${CUSTOM_PCMSolver_LOCATION}/install"
      )
  else()
    ExternalProject_Add(interface_pcmsolver
      PREFIX ${CUSTOM_PCMSolver_LOCATION}
      GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
      GIT_TAG v1.1.3
      CMAKE_ARGS "${PCMSolverCMakeArgs}"
      INSTALL_DIR "${CUSTOM_PCMSolver_LOCATION}/install"
      )
  endif()

  # Set also variables usually set by find_package
  ExternalProject_Get_Property(interface_pcmsolver INSTALL_DIR)
  set(PCMSolver_LIBRARY "${INSTALL_DIR}/lib/libpcm.a")
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include/pcmsolver)  # note [1] below
  set(PCMSolver_INCLUDE_DIRS "${INSTALL_DIR}/include" ${ZLIB_INCLUDE_DIRS})
  set(PCMSolver_LIBRARIES ${PCMSolver_LIBRARY} ${ZLIB_LIBRARIES})

  # Set target for psipcm and psi4 to depend upon as set by find_package
  add_library(PCMSolver::PCMSolver STATIC IMPORTED GLOBAL)
  add_dependencies(PCMSolver::PCMSolver interface_pcmsolver)
  set_target_properties(PCMSolver::PCMSolver PROPERTIES
    IMPORTED_LOCATION "${PCMSolver_LIBRARY}"
    INTERFACE_LINK_LIBRARIES "${PCMSolver_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${PCMSolver_INCLUDE_DIRS}"
    )

  include_directories(SYSTEM "${PCMSolver_INCLUDE_DIRS}")
  set(PCMSolver_PARSE_DIR ${INSTALL_DIR}/bin)
  configure_file(${CMAKE_SOURCE_DIR}/share/python/pcm_placeholder.py.in
    ${CMAKE_SOURCE_DIR}/share/python/pcm_placeholder.py)
endif()

add_definitions(-DHAVE_PCMSOLVER=1)

# [1] It's nice to have a full PCMSolver::PCMSolver target that has embedded
# the library, the linking library paths with dependencies, and the include
# paths with dependencies, just like FindPCMSolver supplies, especially
# since src/bin/libpsipcm needs the headers for compilation and src/bin/psi4
# needs all the linking libraries. Problem is that conventional target
# derived from ExternalProject_Add is a highly sought but not quite
# certified cmake pattern. Hence INTERFACE_INCLUDE_DIRECTORIES complains
# that the directories don't exist at configure time. Hence the hack to
# create an empty directory.
