find_package(xtl 0.6.4 CONFIG QUIET)

if(${xtl_FOUND})
  message(STATUS "${Cyan}Found xtl${ColourReset}: ${xtl_INCLUDE_DIR} (found version ${xtl_VERSION})")
  add_library(xtl_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xtl})
    message(FATAL_ERROR "Suitable xtl could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xtl could not be located, ${Magenta}Building xtl${ColourReset} instead.")
  ExternalProject_Add(xtl_external
    URL
      https://github.com/QuantStack/xtl/archive/0.6.4.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused, but probed by xtl's CMakeLists.txt
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -DXTL_CMAKECONFIG_INSTALL_DIR=share/cmake/xtl # Wrong in upstream CMakeLists.txt
    )

  set(xtl_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xtl CACHE PATH "path to internally built xtlConfig.cmake" FORCE)
endif()
