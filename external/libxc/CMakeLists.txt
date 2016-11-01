# Find/build Libxc

if(FALSE)
#if(${libxc_FOUND})
    message(STATUS "Found libxc: ${libxc_INCLUDE_DIR} (found version ${libxc_VERSION})")
    add_library(libxc_external INTERFACE)  # dummy
else()
  include(ExternalProject)  
  message(STATUS "Suitable libxc could not be located, building one instead.")
  ExternalProject_Add(libxc_external
      URL "http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-3.0.0.tar.gz"
      CONFIGURE_COMMAND "<SOURCE_DIR>/configure"
          "--prefix=${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external"
          "CC=/usr/bin/clang" 
#          "CC=${CMAKE_C_COMPILER}" 
#          "CFLAGS=${CMAKE_C_FLAGS}"
          "--disable-fortran"
          "--bindir=${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/local/bin/"
          "--libdir=${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external/lib/"
          "--includedir=${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/external/include/libxc"
#          "${BUILD_SHARED_LIBS}" 
#          "${ENABLE_STATIC}"
#          "LDFLAGS=${LDFLAGS}"
     BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
     INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )

endif()

