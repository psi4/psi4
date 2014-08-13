# GA_MPICXX_TEST
# --------------
# Attempt to compile a simple MPI program in C++.
#
AC_DEFUN([GA_MPICXX_TEST], [
AS_IF([test "x$with_mpi" != xno], [
    AC_LANG_PUSH([C++])
    AC_CACHE_CHECK([whether a simple C++ MPI program works],
        [ga_cv_cxx_mpi_test], [
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([[#include <mpi.h>]],
[[int myargc; char **myargv; MPI_Init(&myargc, &myargv); MPI_Finalize();]])],
            [ga_cv_cxx_mpi_test=yes],
            [ga_cv_cxx_mpi_test=no])
# That didn't work, so now let's try with our GA_MP_* flags.
        AS_IF([test "x$ga_cv_cxx_mpi_test" = xno], [
        ga_save_LIBS="$LIBS";           LIBS="$LIBS $GA_MP_LIBS"
        ga_save_CPPFLAGS="$CPPFLAGS";   CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
        ga_save_LDFLAGS="$LDFLAGS";     LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([[#include <mpi.h>]],
[[int myargc; char **myargv; MPI_Init(&myargc, &myargv); MPI_Finalize();]])],
            [ga_cv_cxx_mpi_test=yes],
            [ga_cv_cxx_mpi_test=no])
        LIBS="$ga_save_LIBS"
        CPPFLAGS="$ga_save_CPPFLAGS"
        LDFLAGS="$ga_save_LDFLAGS"
        ])
# That didn't work, so now let's try with our GA_MP_* flags and various libs.
        AS_IF([test "x$ga_cv_cxx_mpi_test" = xno], [
        for lib in -lmpi -lmpich; do
        ga_save_LIBS="$LIBS";           LIBS="$LIBS $GA_MP_LIBS $lib"
        ga_save_CPPFLAGS="$CPPFLAGS";   CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
        ga_save_LDFLAGS="$LDFLAGS";     LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([[#include <mpi.h>]],
[[int myargc; char **myargv; MPI_Init(&myargc, &myargv); MPI_Finalize();]])],
            [ga_cv_cxx_mpi_test=$lib; break],
            [ga_cv_cxx_mpi_test=no])
        LIBS="$ga_save_LIBS"
        CPPFLAGS="$ga_save_CPPFLAGS"
        LDFLAGS="$ga_save_LDFLAGS"
        done
        LIBS="$ga_save_LIBS"
        CPPFLAGS="$ga_save_CPPFLAGS"
        LDFLAGS="$ga_save_LDFLAGS"
        ])
        ])
    AC_LANG_POP([C++])
    AS_CASE([$ga_cv_cxx_mpi_test],
        [yes],  [],
        [no],   [AC_MSG_FAILURE([could not link simple C++ MPI program])],
        [*],    [GA_MP_LIBS="$ga_cv_cxx_mpi_test"],
        [])
])
])dnl
