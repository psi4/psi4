# GA_MPIF77_TEST_PROGRAM
# ----------------------
# Create an MPI test program in Fortran 77.
AC_DEFUN([GA_MPIF77_TEST_PROGRAM], [
AC_LANG_PUSH([Fortran 77])
AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
[[      include 'mpif.h'
      integer ierr
      call MPI_Init( ierr )
      call MPI_Finalize( ierr )]])])
AC_LANG_POP([Fortran 77])
])dnl

# GA_MPIF77_TEST_COMPILE
# ----------------------
# Attempt to compile a simple MPI program in Fortran 77.
AC_DEFUN([GA_MPIF77_TEST_COMPILE], [
AC_LANG_PUSH([Fortran 77])
GA_MPIF77_TEST_PROGRAM()
AC_CACHE_CHECK([whether a simple Fortran MPI program compiles],
    [ga_cv_f77_mpi_test_compile],
    [ga_save_FFLAGS="$FFLAGS"; FFLAGS="$FFLAGS $GA_MP_CPPFLAGS"
     AC_COMPILE_IFELSE([],
        [ga_cv_f77_mpi_test_compile=yes],
        [ga_cv_f77_mpi_test_compile=no])
     FFLAGS="$ga_save_FFLAGS"])
rm -f conftest.$ac_ext
AC_LANG_POP([Fortran 77])
AS_IF([test "x$ga_cv_f77_mpi_test_compile" = xno],
    [AC_MSG_FAILURE([could not compile simple Fortran MPI program])])
])dnl

# GA_MPIF77_TEST_LINK
# -------------------
# Attempt to compile a simple MPI program in Fortran 77.
AC_DEFUN([GA_MPIF77_TEST_LINK], [
AC_LANG_PUSH([Fortran 77])
GA_MPIF77_TEST_PROGRAM()
ga_cv_f77_mpi_test_link=no
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([whether a Fortran MPI program links natively])
     AC_LINK_IFELSE([],
        [ga_cv_f77_mpi_test_link=yes
         GA_MP_LIBS=
         GA_MP_LDFLAGS=
         GA_MP_CPPFLAGS=],
        [ga_cv_f77_mpi_test_link=no])
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
# That didn't work.  Let's try adding our GA_MP_* flags.
# The CPPFLAGS are added to FFLAGS since *.f doesn't use CPP.  LIBS changes.
ga_save_LIBS="$LIBS"
ga_save_FFLAGS="$FFLAGS";   FFLAGS="$FFLAGS $GA_MP_CPPFLAGS"
ga_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([whether a Fortran MPI program links with additional env])
     LIBS="$LIBS $GA_MP_LIBS"
     AC_LINK_IFELSE([],
        [ga_cv_f77_mpi_test_link=yes],
        [ga_cv_f77_mpi_test_link=no])
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
# That didn't work, so now let's try with specific libs.
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for mvapich libraries])
     for lib in "-lmpichf90nc -lmpichfarg -lmpich -lpthread" "-lmpichf90 -lmpichfarg -lmpich -pthread" "-lmpichf90nc -lmpichfarg -lmpich" "-lmpichf90 -lmpichfarg -lmpich" "-lmpichfarg -lmpich -lpthread" "-lmpichfarg -lmpich"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_f77_mpi_test_link="$lib"; break],
            [ga_cv_f77_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for mpich libraries])
     for lib in "-lmpichf90 -lmpich -lpthread" "-lmpichf90 -lmpich" "-lmpich -pthread" "-lmpich"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_f77_mpi_test_link="$lib"; break],
            [ga_cv_f77_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for hpmpi libraries])
     for lib in "-lhpmpio -lhpmpi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_f77_mpi_test_link="$lib"; break],
            [ga_cv_f77_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for intelmpi libraries])
     for lib in "-lmpi -lmpigf -lmpigi -lpthread" "-lmpi -lmpigf -lmpigi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_f77_mpi_test_link="$lib"; break],
            [ga_cv_f77_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
AS_IF([test "x$ga_cv_f77_mpi_test_link" = xno],
    [AC_MSG_CHECKING([for openmpi libraries])
     for lib in "-lmpi_f90 -lmpi_f77 -lmpi" "-lmpi_f77 -lmpi"
     do
        LIBS="$LIBS $lib"
        AC_LINK_IFELSE([],
            [ga_cv_f77_mpi_test_link="$lib"; break],
            [ga_cv_f77_mpi_test_link=no])
        LIBS="$ga_save_LIBS"
     done
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_cv_f77_mpi_test_link])])
rm -f conftest.$ac_ext
LIBS="$ga_save_LIBS"
LDFLAGS="$ga_save_LDFLAGS"
FFLAGS="$ga_save_FFLAGS"
AC_LANG_POP([Fortran 77])
AS_CASE([$ga_cv_f77_mpi_test_link],
    [yes],  [],
    [no],   [AC_MSG_FAILURE([could not link simple Fortran MPI program])],
    [*],    [GA_MP_LIBS="$ga_cv_f77_mpi_test_link"],
    [])
])dnl
