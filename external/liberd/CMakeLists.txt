include(ExternalProject)

find_package(liberd QUIET)
if((NOT liberd_FOUND) AND ENABLE_LIBERD)
   if(BUILD_FPIC)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
   endif()
   message(STATUS 
      "Suitable Liberd could not be located, building one instead")
   ExternalProject_Add(liberd_external
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/liberd-src
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/external/liberd
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                 -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
     INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install 
                     DESTDIR=${CMAKE_BINARY_DIR}/stage
      ) 
else()
    add_library(liberd_external INTERFACE)
endif()

