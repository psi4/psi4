# GA_MPICC_TEST_PROGRAM
# ---------------------
# Create an MPI test program in C.
AC_DEFUN([GA_MPICC_TEST_PROGRAM], [
AC_LANG_PUSH([C])
AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <mpi.h>]],
[[int myargc; char **myargv; MPI_Init(&myargc, &myargv); MPI_Finalize();]])])
AC_LANG_POP([C])
])dnl

# GA_MPICC_TEST_COMPILE
# ---------------------
# Attempt to compile a simple MPI program in C.
AC_DEFUN([GA_MPICC_TEST_COMPILE], [
AC_LANG_PUSH([C])
GA_MPICC_TEST_PROGRAM()
AC_CACHE_CHECK([whether a simple C MPI program compiles],
    [ga_cv_c_mpi_test_compile],
    [ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
     AC_COMPILE_IFELSE([],
        [ga_cv_c_mpi_test_compile=yes],
        [ga_cv_c_mpi_test_compile=no])
     CPPFLAGS="$ga_save_CPPFLAGS"])
rm -f conftest.$ac_ext
AC_LANG_POP([C])
AS_IF([test "x$ga_cv_c_mpi_test_compile" = xno],
    [AC_MSG_FAILURE([could not compile simple C MPI program])])
])dnl

# GA_MPICC_TEST_LINK
# ------------------
# Attempt to link a simple MPI program in C.
AC_DEFUN([GA_MPICC_TEST_LINK], [
AC_LANG_PUSH([C])
GA_MPICC_TEST_PROGRAM()
ga_cv_c_mpi_test_link=no
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([whether a C MPI program links natively])
     AC_LINK_IFELSE([],
        [ga_cv_c_mpi_test_link=yes
         GA_MP_LIBS=
         GA_MP_LDFLAGS=
         GA_MP_CPPFLAGS=],
        [ga_cv_c_mpi_test_link=no])
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
# That didn't work, so now let's try adding our GA_MP_* flags.
# The CPPFLAGS and LDFLAGS are added up top here, but LIBS will change.
ga_save_LIBS="$LIBS"
ga_save_CPPFLAGS="$CPPFLAGS";  CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
ga_save_LDFLAGS="$LDFLAGS";    LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([whether a C MPI program links with additional env])
     LIBS="$LIBS $GA_MP_LIBS"
     AC_LINK_IFELSE([],
        [ga_cv_c_mpi_test_link=yes],
        [ga_cv_c_mpi_test_link=no])
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
# That didn't work, so now let's try with specific libs.
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for mvapich libraries])
     for lib in "-lmpich -lpthread" "-lmpich"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_c_mpi_test_link="$lib"; break],
            [ga_cv_c_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for mpich libraries])
     for lib in "-lmpich -lpthread" "-lmpich"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_c_mpi_test_link="$lib"; break],
            [ga_cv_c_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for hpmpi libraries])
     for lib in "-lhpmpio -lhpmpi" "-lhpmpi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_c_mpi_test_link="$lib"; break],
            [ga_cv_c_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for intelmpi libraries])
     for lib in "-lmpi -lmpiif -lmpigi -lrt -lpthread" "-lmpi -lmpiif -lmpigi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_c_mpi_test_link="$lib"; break],
            [ga_cv_c_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
AS_IF([test "x$ga_cv_c_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for openmpi libraries])
     for lib in "-lmpi -lpthread" "-lmpi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_c_mpi_test_link="$lib"; break],
            [ga_cv_c_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_c_mpi_test_link])])
rm -f conftest.$ac_ext
LIBS="$ga_save_LIBS"
LDFLAGS="$ga_save_LDFLAGS"
CPPFLAGS="$ga_save_CPPFLAGS"
AC_LANG_POP([C])
AS_CASE([$ga_cv_c_mpi_test_link],
    [yes],  [],
    [no],   [AC_MSG_FAILURE([could not link a C MPI program])],
    [*],    [GA_MP_LIBS="$ga_cv_c_mpi_test_link"],
    [])
])dnl
