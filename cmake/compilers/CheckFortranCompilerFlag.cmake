#.rst:
# CheckCCompilerFlag
# ------------------
#
# Check whether the C compiler supports a given flag.
#
# CHECK_C_COMPILER_FLAG(<flag> <var>)
#
# ::
#
#   <flag> - the compiler flag
#   <var>  - variable to store the result
#            Will be created as an internal cache variable.
#
# This internally calls the check_c_source_compiles macro and sets
# CMAKE_REQUIRED_DEFINITIONS to <flag>.  See help for
# CheckCSourceCompiles for a listing of variables that can otherwise
# modify the build.  The result only tells that the compiler does not
# give an error message when it encounters the flag.  If the flag has
# any effect or even a specific one is beyond the scope of this module.

#=============================================================================
# Copyright 2006-2011 Kitware, Inc.
# Copyright 2006 Alexander Neundorf <neundorf@kde.org>
# Copyright 2011 Matthias Kretz <kretz@kde.org>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

include(CheckFortranSourceCompiles)
include(CMakeCheckCompilerFlagCommonPatterns)

macro (CHECK_Fortran_COMPILER_FLAG _FLAG _RESULT)
   set(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
   set(CMAKE_REQUIRED_DEFINITIONS "${_FLAG}")

   # Normalize locale during test compilation.
   set(_CheckFortranCompilerFlag_LOCALE_VARS LC_ALL LC_MESSAGES LANG)
   foreach(v ${_CheckFortranCompilerFlag_LOCALE_VARS})
     set(_CheckFortranCompilerFlag_SAVED_${v} "$ENV{${v}}")
     set(ENV{${v}} Fortran)
   endforeach()
   CHECK_COMPILER_FLAG_COMMON_PATTERNS(_CheckFortranCompilerFlag_COMMON_PATTERNS)
   # Don't f**k up the indentation here!!! CMake expects fixed format Fortran sources!!!
   set(TEST_SOURCE "       program test\n       end program test") 
   CHECK_Fortran_SOURCE_COMPILES("${TEST_SOURCE}" ${_RESULT}
     # Some compilers do not fail with a bad flag
     FAIL_REGEX "command line option .* is valid for .* but not for C" # GNU
     ${_CheckFortranCompilerFlag_COMMON_PATTERNS}
     )
   foreach(v ${_CheckFortranCompilerFlag_LOCALE_VARS})
     set(ENV{${v}} ${_CheckFortranCompilerFlag_SAVED_${v}})
     unset(_CheckFortranCompilerFlag_SAVED_${v})
   endforeach()
   unset(_CheckFortranCompilerFlag_LOCALE_VARS)
   unset(_CheckFortranCompilerFlag_COMMON_PATTERNS)

   set (CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
endmacro ()
