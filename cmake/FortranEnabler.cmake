# This macro enables Fortran as a language for the Psi project.
# It relies on the FortranCInterface CMake module and will:
# 1. check that C/C++ and Fortran ABIs are compatible. CMake will issue
#    a fatal error if they are not;
# 2. generate the FCMangle.h header. This header contains macros
#    to perform the correct name mangling of Fortran/C interfacing
#    subroutines/functions.
# The macro might be called in different places, at different times.
# This might cause multiple, redundant, calls to the FortranCInterface_VERIFY(CXX) function
# and multiple, redundant, generation of the FCMangle.h header.
# To avoid this, the cache variable FORTRAN_ENABLER_CALLED is set once the macro is called
# for the first time. 
# The macro checks if the variable is set to TRUE. If yes, execution is skipped.

macro(fortran_enabler)
    if(NOT FORTRAN_ENABLER_CALLED)
       if(CMAKE_Fortran_COMPILER)
          # Determine Fortran name mangling, used for external linking                             
          if ((${CMAKE_SYSTEM_NAME} MATCHES "Darwin") AND (CMAKE_SYSTEM_VERSION VERSION_LESS 13))
              # 12.5.0 broken (LAB), above should trap pre-Mavericks
              set(CMAKE_C_FLAGS_HOLD ${CMAKE_C_FLAGS})
              set(CMAKE_CXX_FLAGS_HOLD ${CMAKE_CXX_FLAGS})
              message(STATUS "Trying Fortran name mangling appending -m32 flag")
              set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
              set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
          endif()
                                                                                                  
          include(FortranCInterface)
          FortranCInterface_VERIFY(CXX)
          init_FCMangle()
                                                                                                  
          if ((${CMAKE_SYSTEM_NAME} MATCHES "Darwin") AND (CMAKE_SYSTEM_VERSION VERSION_LESS 13))
              set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_HOLD})
              set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_HOLD})
          endif()
          
	  # Set FORTRAN_ENABLER_CALLED to TRUE
	  set(FORTRAN_ENABLER_CALLED TRUE CACHE BOOL "Fortran enabler called") 
       else(CMAKE_Fortran_COMPILER)
            message(FATAL_ERROR "Compilation of Fortran sources requires a Fortran compiler!")
       endif(CMAKE_Fortran_COMPILER)
    endif()
endmacro(fortran_enabler)
