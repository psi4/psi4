include(ExternalProject)
option(ENABLE_BOUDNS_CHECK "Enables bounds check in Fortran" OFF)
find_package(dkh QUIET)

if(NOT DKH_FOUND AND ${ENABLE_DKH})
   if(BUILD_FPIC)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
   endif()
   message(STATUS "DKH not found. The pre-packaged version will be build.")
   ExternalProject_Add(dkh_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dkh-src
      CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                 -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/dkh
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                      DESTDIR=${CMAKE_BINARY_DIR}/stage
   )
else()
   add_library(dkh_external INTERFACE)
endif()

