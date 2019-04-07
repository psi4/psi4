find_package(xtensor-blas 0.15.2 CONFIG QUIET)

if(${xtensor-blas_FOUND})
  message(STATUS "${Cyan}Found xtensor-blas${ColourReset}: ${xtensor-blas_INCLUDE_DIR} (found version ${xtensor-blas_VERSION})")
  add_library(xtensor-blas_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xtensor-blas})
    message(FATAL_ERROR "Suitable xtensor-blas could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xtensor-blas could not be located, ${Magenta}Building xtensor-blas${ColourReset} instead.")
  ExternalProject_Add(xtensor-blas_external
    DEPENDS
      xtensor_external
    URL
      https://github.com/QuantStack/xtensor-blas/archive/0.16.0.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused, but probed by xtensor-blas's CMakeLists.txt
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -Dxtensor_DIR=${xtensor_DIR}
      -DXTENSOR_USE_FLENS_BLAS=OFF # Turn ON if no suitable BLAS/LAPACK found?
      -DXTENSOR_BLAS_CMAKECONFIG_INSTALL_DIR=share/cmake/xtensor-blas # Wrong in upstream CMakeLists.txt
    )

  set(xtensor-blas_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xtensor-blas CACHE PATH "path to internally built xtensor-blasConfig.cmake" FORCE)
endif()
