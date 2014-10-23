# Based on the VcMacros.cmake file from the Vc library. Adapted and significantly
# simplified to get it to work in PCMSolver (since it uses a different compiler
# setup strategy)
# Macros for use with the Vc library. Vc can be found at http://code.compeng.uni-frankfurt.de/projects/vc
#
# The following macros are provided:
# vc_determine_compiler
# vc_set_preferred_compiler_flags
#
#=============================================================================
# Copyright 2009-2013   Matthias Kretz <kretz@kde.org>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * The names of Kitware, Inc., the Insight Consortium, or the names of
#    any consortium members, or of any contributors, may not be used to
#    endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# Based on vc_check_assembler from the Vc project
macro(check_assembler)
   if(APPLE)
      if(NOT CMAKE_CXX_COMPILER_ID MATCHES Clang)
         message(WARNING "Apple does not provide an assembler with AVX support. AVX will not be available. Please use Clang if you want to use AVX.")
         set(DEFINITIONS "${DEFINITIONS} -DNO_XGETBV")
         set(AVX_INTRINSICS_BROKEN true)
      endif()
   else(APPLE)
      if(${ARGC} EQUAL 1)
         set(_as "${ARGV1}")
      else()
         exec_program(${CMAKE_CXX_COMPILER} ARGS -print-prog-name=as OUTPUT_VARIABLE _as)
         mark_as_advanced(_as)
      endif()
      if(NOT _as)
         message(WARNING "Could not find 'as', the assembler used by GCC. Hoping everything will work out...")
      else()
         exec_program(${_as} ARGS --version OUTPUT_VARIABLE _as_version)
         string(REGEX REPLACE "\\([^\\)]*\\)" "" _as_version "${_as_version}")
         string(REGEX MATCH "[1-9]\\.[0-9]+(\\.[0-9]+)?" _as_version "${_as_version}")
         if(_as_version VERSION_LESS "2.18.93")
	    message(WARNING "Your binutils is too old (${_as_version}). Some optimizations of will be disabled.")
            add_definitions(-DNO_XGETBV) # old assembler doesn't know the xgetbv instruction
         endif()
      endif()
   endif(APPLE)
endmacro()

# Based on vc_check_fpmath
macro(check_fpmath)
   # if compiling for 32 bit x86 we need to use the -mfpmath=sse since the x87 is broken by design
   include (CheckCXXSourceRuns)
   check_cxx_source_runs("int main() { return sizeof(void*) != 8; }" VOID_PTR_IS_64BIT)
   if(NOT VOID_PTR_IS_64BIT)
      exec_program(${CMAKE_C_COMPILER} ARGS -dumpmachine OUTPUT_VARIABLE _gcc_machine)
      if(_gcc_machine MATCHES "[x34567]86" OR _gcc_machine STREQUAL "mingw32")
	 set(DEFINTIONS "${DEFINITIONS} -mfpmath=sse")     
	      #vc_add_compiler_flag(DEFINITIONS "-mfpmath=sse")
      endif()
   endif()
endmacro()

# Based on vc_set_preferred_compiler_flags
macro(find_broken_intrinsics)

   set(SSE_INTRINSICS_BROKEN false)
   set(AVX_INTRINSICS_BROKEN false)
   set(XOP_INTRINSICS_BROKEN false)
   set(FMA4_INTRINSICS_BROKEN false)

   if(CMAKE_CXX_COMPILER_ID MATCHES GCC)
      set(DEFINITIONS "${DEFINITIONS} -Wabi -fabi-version=0") # ABI version 4 is required to make __m128 and __m256 appear as different types. 0 should give us the latest version.

      # GCC 4.5.[01] fail at inlining some functions, creating functions with a single instructions,
      # thus creating a large overhead.
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.5.2" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.5.0")
	 message(WARNING "GCC 4.5.0 and 4.5.1 have problems with inlining correctly. Setting early-inlining-insns=12 as workaround.")
	 set(DEFINITIONS "${DEFINITIONS} --param early-inlining-insns=12")
      endif()

      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.1.99")
	 message(WARNING "Your GCC is ancient and crashes on some important optimizations.\n The full set of SSE2 intrinsics is not supported.\n PCMSolver will fall back to the scalar implementation.\n Use of the may_alias and always_inline attributes will be disabled.\n In turn all code using PCMSolver must be compiled with -fno-strict-aliasing")
	 set(DEFINTIONS "${DEFINITIONS} -fno-strict-aliasing")
         set(AVX_INTRINSICS_BROKEN true)
         set(SSE_INTRINSICS_BROKEN true)
      elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.4.6")
	 message(WARNING "Your GCC is older than 4.4.6. This is known to cause problems/bugs. Please update to the latest GCC if you can.")
         set(AVX_INTRINSICS_BROKEN true)
	 if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.3.0")
            message(WARNING "Your GCC is older than 4.3.0. It is unable to handle the full set of SSE2 intrinsics. All SSE code will be disabled. Please update to the latest GCC if you can.")
            set(SSE_INTRINSICS_BROKEN true)
         endif()
      endif()

      if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.6.0)
         message(WARNING "GCC 4.6.0 miscompiles AVX loads/stores, leading to spurious segfaults. Disabling AVX per default.")
         set(AVX_INTRINSICS_BROKEN true)
      elseif(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.7.0)
         message(WARNING "GCC 4.7.0 miscompiles at -O3, adding -fno-predictive-commoning to the compiler flags as workaround")
         set(DEFINITIONS "${DEFINITIONS} -fno-predictive-commoning")
      elseif(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.8.0)
         message(WARNING "GCC 4.8.0 miscompiles at -O3, adding -fdisable-tree-vrp2 to the compiler flags as workaround")
         set(DEFINITIONS "${DEFINITIONS} -fdisable-tree-vrp2")
      endif()

      check_fpmath()
      check_assembler()
   elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
      # Intel doesn't implement the XOP or FMA4 intrinsics
      set(XOP_INTRINSICS_BROKEN true)
      set(FMA4_INTRINSICS_BROKEN true)
   elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
         # the LLVM assembler gets FMAs wrong (bug 15040)
	 set(DEFINITIONS "${DEFINITIONS} -no-integrated-as")
      endif()

   endif()

endmacro()
