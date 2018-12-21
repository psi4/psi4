find_package(xsimd 7.1.3 CONFIG QUIET)

if(${xsimd_FOUND})
  message(STATUS "${Cyan}Found xsimd${ColourReset}: ${xsimd_INCLUDE_DIR} (found version ${xsimd_VERSION})")
  add_library(xsimd_external INTERFACE)  # dummy
else()
  if(${CMAKE_INSIST_FIND_PACKAGE_xsimd})
    message(FATAL_ERROR "Suitable xsimd could not be externally located as user insists")
  endif()

  include(ExternalProject)
  message(STATUS "Suitable xsimd could not be located, ${Magenta}Building xsimd${ColourReset} instead.")
  ExternalProject_Add(xsimd_external
    URL
      https://github.com/QuantStack/xsimd/archive/7.1.3.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused, but probed by xsimd's CMakeLists.txt
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
      -DCMAKE_CXX_STANDARD=${psi4_CXX_STANDARD}
      -DXSIMD_CMAKECONFIG_INSTALL_DIR=share/cmake/xsimd # Wrong in upstream CMakeLists.txt
    )

  set(xsimd_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/xsimd CACHE PATH "path to internally built xsimdConfig.cmake" FORCE)
endif()
