find_package(xtensor 0.18.3 CONFIG QUIET)

if(${xtensor_FOUND})
  message(STATUS "${Cyan}Found xtensor${ColourReset}: ${xtensor_INCLUDE_DIR} (found version ${xtensor_VERSION})")
  add_library(xtensor_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xtensor})
    message(FATAL_ERROR "Suitable xtensor could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xtensor could not be located, ${Magenta}Building xtensor${ColourReset} instead.")
  ExternalProject_Add(xtensor_external
    DEPENDS
      xtl_external
      xsimd_external
    GIT_REPOSITORY
      https://github.com/QuantStack/xtensor
    GIT_TAG
      0.18.3
    UPDATE_COMMAND
      ""
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused, but probed by xtensor's CMakeLists.txt
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -Dxtl_DIR=${xtl_DIR}
      -Dxsimd_DIR=${xsimd_DIR}
      -DXTENSOR_USE_XSIMD=ON
      -DDOWNLOAD_GBENCHMARK=OFF
      -DXTENSOR_CMAKECONFIG_INSTALL_DIR=share/cmake/xtensor # Wrong in upstream CMakeLists.txt
    )

  set(xtensor_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xtensor CACHE PATH "path to internally built xtensorConfig.cmake" FORCE)
endif()
