find_package(HighFive ${HighFive_pinned} CONFIG QUIET)
if(${HighFive_FOUND})
  message(STATUS "${Cyan}Found HighFive${ColourReset}: ${HighFive_INCLUDE_DIR} (found version ${HighFive_VERSION})")
  add_library(HighFive_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_HighFive})
    message(FATAL_ERROR "Suitable HighFive could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable HighFive could not be located, ${Magenta}Building HighFive${ColourReset} instead.")
  ExternalProject_Add(HighFive_external
    URL
      https://github.com/BlueBrain/HighFive/archive/v2.2.1.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -DHIGHFIVE_USE_BOOST=OFF
      -DHIGHFIVE_USE_EIGEN=OFF
      -DHIGHFIVE_EXAMPLES=OFF
      -DHIGHFIVE_UNIT_TESTS=OFF
      -DHIGHFIVE_USE_INSTALL_DEPS=ON
    )

  set(HighFive_DIR ${STAGED_INSTALL_PREFIX}/share/HighFive/CMake CACHE PATH "path to internally built HighFiveConfig.cmake" FORCE)
endif()
