#!/bin/sh

CCOMP=${HOME}/bin/mpicc
CXXCOMP=${HOME}/bin/mpicxx
FORTCOMP=${HOME}/bin/mpif77
BOOSTROOT=${HOME}
BLASLIB=${HOME}/lib/mkl
BLASINC=${HOME}/include/mkl
LAPACKLIB=${BLASLIB}
LAPACKINC=${BLASINC}
SCALAPACKLIB=${BLASLIB}
SCALAPACKINC=${BLASINC}
BUILDDIR="build"
#GAINC=include
#GALIB=lib
if [ ! -d "${BUILDDIR}" ];then
   mkdir "${BUILDDIR}"
fi

cd ${BUILDDIR}
../configure.cmake --with-cc=${CCOMP} --with-cxx=${CXXCOMP}\
            --with-f77=${FORTCOMP} --with-boost-root=${BOOSTROOT} \
            --with-blas-lib="-L${BLASLIB} -lmkl" --with-blas-inc=${BLASINC} \
            --with-verbose
 #           --with-ga-inc=${GAINC} --with-ga-lib=${GALIB} \
make
