#Try to find MKL first
message(STATUS ${OSX_ARCHITECTURES})
find_path(BLAS_INCLUDE_DIR mkl_cblas.h
          PATHS $ENV{MKLROOT}/include
)
find_path(LAPACK_INCLUDE_DIR mkl_lapacke.h
          PATHS $ENV{MKLROOT}/include
)
find_library(BLAS_LIBRARY mkl_core
           PATHS
)


