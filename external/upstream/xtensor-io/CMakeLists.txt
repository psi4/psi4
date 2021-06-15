find_package(xtensor-io ${xtensor-io_pinned} CONFIG QUIET)
if(${xtensor-io_FOUND})
  message(STATUS "${Cyan}Found xtensor-io${ColourReset}: ${xtensor-io_INCLUDE_DIR} (found version ${xtensor-io_VERSION})")
  add_library(xtensor-io_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xtensor-io})
    message(FATAL_ERROR "Suitable xtensor-io could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xtensor-io could not be located, ${Magenta}Building xtensor-io${ColourReset} instead.")
  ExternalProject_Add(xtensor-io_external
    DEPENDS
      xtl_external
      xtensor_external
      HighFive_external
    URL
      https://github.com/QuantStack/xtensor-io/archive/${xtensor-io_pinned}.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -DDOWNLOAD_GBENCHMARK=OFF
      -Dxtl_DIR=${xtl_DIR}
      -Dxtensor_DIR=${xtensor_DIR}
      -DHighFive_DIR=${HighFive_DIR}
      -DXTENSOR_IO_CMAKECONFIG_INSTALL_DIR=share/cmake/xtensor-io # Wrong in upstream CMakeLists.txt
    PATCH_COMMAND
      patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/highfive-dumpmode.patch
    )

  set(xtensor-io_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xtensor-io CACHE PATH "path to internally built xtensor-ioConfig.cmake" FORCE)
endif()
