include(ExternalProject)
option(ENABLE_BOUDNS_CHECK "Enables bounds check in Fortran" OFF)

if(NOT DKH_FOUND AND ENABLE_DKH)
   ExternalProject_Add(dkh_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dkh-src
      CMAKE_ARGS -CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                 -CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                      DESTDIR=${CMAKE_BINARY_DIR}/stage
   )
else()
   add_library(dkh_external INTERFACE)
endif()

