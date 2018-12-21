find_package(xtensor-python 0.22.1 CONFIG QUIET)

if(${xtensor-python_FOUND})
  message(STATUS "${Cyan}Found xtensor-python${ColourReset}: ${xtensor-python_INCLUDE_DIR} (found version ${xtensor-python_VERSION})")
  add_library(xtensor-python_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xtensor-python})
    message(FATAL_ERROR "Suitable xtensor-python could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xtensor-python could not be located, ${Magenta}Building xtensor-python${ColourReset} instead.")
  ExternalProject_Add(xtensor-python_external
    DEPENDS
      pybind11_external
      xtl_external
      xtensor_external
    URL
      https://github.com/QuantStack/xtensor-python/archive/0.22.1.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused, but probed by xtensor-python's CMakeLists.txt
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -Dpybind11_DIR=${pybind11_DIR}
      -Dxtl_DIR=${xtl_DIR}
      -Dxtensor_DIR=${xtensor_DIR}
      -DXTENSOR_PYTHON_CMAKECONFIG_INSTALL_DIR=share/cmake/xtensor-python # Wrong in upstream CMakeLists.txt
    )

  set(xtensor-python_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xtensor-python CACHE PATH "path to internally built xtensor-pythonConfig.cmake" FORCE)
endif()
