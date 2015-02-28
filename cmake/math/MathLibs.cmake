

# -----------------------------------------------------------------------------
# Copyright 2011-2013 Jonas Juselius <firstname.lastname at uit.no>
#                     Radovan Bast   <lastname at kth.se>
#
# Distributed under the GNU Lesser General Public License.
# See accompanying file LICENSE for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# -----------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SYSTEM_NATIVE

set(SYSTEM_NATIVE_BLAS_INCLUDE_PATH_SUFFIXES)
set(SYSTEM_NATIVE_LAPACK_INCLUDE_PATH_SUFFIXES)

set(SYSTEM_NATIVE_BLAS_HEADERS   cblas.h)
set(SYSTEM_NATIVE_LAPACK_HEADERS clapack.h)

set(SYSTEM_NATIVE_BLAS_LIBRARY_PATH_SUFFIXES)
set(SYSTEM_NATIVE_LAPACK_LIBRARY_PATH_SUFFIXES)

set(SYSTEM_NATIVE_BLAS_LIBS   blas)
set(SYSTEM_NATIVE_LAPACK_LIBS lapack)

#-------------------------------------------------------------------------------
# ESSL

set(ESSL_BLAS_INCLUDE_PATH_SUFFIXES)
set(ESSL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ESSL_BLAS_HEADERS   UNKNOWN)
set(ESSL_LAPACK_HEADERS UNKNOWN)

set(ESSL_BLAS_LIBRARY_PATH_SUFFIXES)
set(ESSL_LAPACK_LIBRARY_PATH_SUFFIXES)

set(ESSL_BLAS_LIBS   essl)
set(ESSL_LAPACK_LIBS essl)

#-------------------------------------------------------------------------------
# ACML

set(ACML_BLAS_INCLUDE_PATH_SUFFIXES)
set(ACML_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ACML_BLAS_HEADERS   cblas.h)
set(ACML_LAPACK_HEADERS clapack.h)

set(ACML_BLAS_LIBRARY_PATH_SUFFIXES   libso)
set(ACML_LAPACK_LIBRARY_PATH_SUFFIXES libso)

set(ACML_BLAS_LIBS   acml)
set(ACML_LAPACK_LIBS acml)

#-------------------------------------------------------------------------------
# ATLAS

set(ATLAS_BLAS_INCLUDE_PATH_SUFFIXES   atlas)
set(ATLAS_LAPACK_INCLUDE_PATH_SUFFIXES atlas)

set(ATLAS_BLAS_HEADERS   cblas.h)
set(ATLAS_LAPACK_HEADERS clapack.h)

set(ATLAS_BLAS_LIBRARY_PATH_SUFFIXES   atlas atlas-base atlas-base/atlas atlas-sse3)
set(ATLAS_LAPACK_LIBRARY_PATH_SUFFIXES atlas atlas-base atlas-base/atlas atlas-sse3)

set(ATLAS_BLAS_LIBS   f77blas cblas atlas)
set(ATLAS_LAPACK_LIBS atlas lapack)

#-------------------------------------------------------------------------------
# MKL

set(MKL_BLAS_INCLUDE_PATH_SUFFIXES)
set(MKL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(MKL_BLAS_HEADERS   mkl_cblas.h)
set(MKL_LAPACK_HEADERS mkl_lapack.h)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   intel64 em64t)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES intel64 em64t)
else()
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   ia32 32)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES ia32 32)
endif()

set(_thread_lib)
if(ENABLE_THREADED_MKL)
    if(MKL_COMPILER_BINDINGS MATCHES Intel)
        set(_thread_lib mkl_intel_thread)
    endif()
    if(MKL_COMPILER_BINDINGS MATCHES PGI)
        set(_thread_lib mkl_pgi_thread)
    endif()
    if(MKL_COMPILER_BINDINGS MATCHES GNU)
        set(_thread_lib mkl_gnu_thread)
    endif()
    if(MKL_COMPILER_BINDINGS MATCHES Clang)
        set(_thread_lib mkl_gnu_thread)
    endif()
else()
    set(_thread_lib mkl_sequential)
endif()

if(MKL_COMPILER_BINDINGS MATCHES Intel)
    set(_compiler_mkl_interface mkl_intel)
endif()
if(MKL_COMPILER_BINDINGS MATCHES PGI)
    set(_compiler_mkl_interface mkl_intel)
endif()
if(MKL_COMPILER_BINDINGS MATCHES GNU)
    set(_compiler_mkl_interface mkl_gf)
endif()
if(MKL_COMPILER_BINDINGS MATCHES Clang)
    set(_compiler_mkl_interface mkl_gf)
endif()

set(_lib_suffix)
if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    if(ENABLE_64BIT_INTEGERS)
        set(_lib_suffix _ilp64)
    else()
        set(_lib_suffix _lp64)
    endif()
endif()

if(ENABLE_SCALAPACK)
    set(_scalapack_lib      mkl_scalapack${_lib_suffix})
    if(${BLACS_IMPLEMENTATION} STREQUAL "intelmpi")
        set(_blacs_lib mkl_blacs_intelmpi${_lib_suffix})
    elseif(${BLACS_IMPLEMENTATION} STREQUAL "openmpi")
        set(_blacs_lib mkl_blacs_openmpi${_lib_suffix})
    elseif(${BLACS_IMPLEMENTATION} STREQUAL "sgimpt")
        set(_blacs_lib mkl_blacs_sgimpt${_lib_suffix})
    else()
        message(FATAL_ERROR "BLACS implementation ${BLACS_IMPLEMENTATION} not recognized/supported")
    endif()
else()
    set(_scalapack_lib)
    set(_blacs_lib)
endif()

# miro: for MKL 10.0.1.014
set(MKL_BLAS_LIBS ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core mkl_def mkl_mc ${_blacs_lib} guide pthread m)
#  try this MKL BLAS combination with SGI MPT
set(MKL_BLAS_LIBS2 ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_lib}   guide pthread m)
# newer MKL BLAS versions do not have libguide
set(MKL_BLAS_LIBS3 ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_lib}         pthread m)
# ancient MKL BLAS
set(MKL_BLAS_LIBS4 mkl guide m)

set(MKL_LAPACK_LIBS mkl_lapack95${_lib_suffix} ${_compiler_mkl_interface}${_lib_suffix})

# older MKL LAPACK
set(MKL_LAPACK_LIBS2 mkl_lapack)

unset(_lib_suffix)
unset(_thread_lib)
unset(_compiler_mkl_interface)
unset(_scalapack_lib)
unset(_blacs_lib)
