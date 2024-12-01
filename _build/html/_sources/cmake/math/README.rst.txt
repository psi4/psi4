

Where does CMake search math libraries if you specify --blas/lapack=auto?
-------------------------------------------------------------------------

CMake will look in the environment variable MATH_ROOT.

For instance my .bashrc contains::

  source /opt/intel/bin/compilervars.sh intel64
  export MATH_ROOT=/opt/intel/mkl


Order of math libraries
-----------------------

Order is set by MATH_LIB_SEARCH_ORDER in MathLibs.cmake.
You can override this order by setting BLAS_TYPE and/or LAPACK_TYPE
for example to ATLAS or some other library that you prefer.


What to edit if your math library is not found although you have set MATH_ROOT?
-------------------------------------------------------------------------------

Normally you only need to edit MathLibs.cmake to add new libraries
or edit existing ones.

Since a vendor can provide libraries with different "fingerprints"
(example MKL), you can define different combinations (up to 9), for instance::

  set(MKL_BLAS_LIBS  ...)
  set(MKL_BLAS_LIBS2 ...)
  set(MKL_BLAS_LIBS3 ...)
  set(MKL_BLAS_LIBS4 ...)
  set(MKL_BLAS_LIBS5 ...)

Then CMake will first try MKL_BLAS_LIBS, then MKL_BLAS_LIBS2, etc.
The first pattern that will match will be linked against.
