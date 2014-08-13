# GA_DISABLE_MPI_TESTS()
# ---------------------
# Whether to disable all MPI linker tests.
AC_DEFUN([GA_DISABLE_MPI_TESTS],
[AC_ARG_ENABLE([mpi-tests],
    [AS_HELP_STRING([--disable-mpi-tests], [disable MPI linker tests])],
    [],
    [enable_mpi_tests=yes])
])# GA_DISABLE_MPI_TESTS
