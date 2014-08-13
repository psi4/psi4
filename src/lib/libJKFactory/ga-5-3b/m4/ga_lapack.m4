# GA_F77_LAPACK_TEST
# ------------------
# Generate Fortran 77 conftest for LAPACK.
AC_DEFUN([GA_F77_LAPACK_TEST], [AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
[[      implicit none
      external DGETRS
      CALL DGETRS ()]])])
])


# GA_C_LAPACK_TEST
# ----------------
# Generate C conftest for LAPACK.
AC_DEFUN([GA_C_LAPACK_TEST], [AC_LANG_CONFTEST([AC_LANG_PROGRAM(
[#ifdef __cplusplus
extern "C" {
#endif
char dgetrs ();
#ifdef __cplusplus
}
#endif
],
[[char result = dgetrs ();
]])])
])


# GA_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# ---------------------------------------------------
# Test for LAPACK. Modelled after GA_BLAS macro.
AC_DEFUN([GA_LAPACK],
[AC_REQUIRE([GA_BLAS])
AC_ARG_WITH([lapack],
    [AS_HELP_STRING([--with-lapack=[[ARG]]], [use external LAPACK library])])

ga_lapack_ok=no
AS_IF([test "x$with_lapack" = xno], [ga_lapack_ok=skip])

# Parse --with-lapack argument. Clear previous values first.
LAPACK_LIBS=
LAPACK_LDFLAGS=
LAPACK_CPPFLAGS=
GA_ARG_PARSE([with_lapack], [LAPACK_LIBS], [LAPACK_LDFLAGS], [LAPACK_CPPFLAGS])

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(dgetrs)

ga_save_LIBS="$LIBS"
ga_save_LDFLAGS="$LDFLAGS"
ga_save_CPPFLAGS="$CPPFLAGS"

LDFLAGS="$LAPACK_LDFLAGS $BLAS_LDFLAGS $LDFLAGS"
CPPFLAGS="$LAPACK_CPPFLAGS $BLAS_CPPFLAGS $CPPFLAGS"

AC_MSG_NOTICE([Attempting to locate LAPACK library])

# We cannot use LAPACK if BLAS is not found
AS_IF([test $ga_blas_ok != yes], [ga_lapack_ok=noblas])

# First, check environment/command-line variables.
# If failed, erase LAPACK_LIBS but maintain LAPACK_LDFLAGS and LAPACK_CPPFLAGS.
AS_IF([test $ga_lapack_ok = no],
    [LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
     AS_IF([test "x$enable_f77" = xno],
        [AC_MSG_CHECKING([for C LAPACK with user-supplied flags])
         AC_LANG_PUSH([C])
         GA_C_LAPACK_TEST()
         AC_LINK_IFELSE([], [ga_lapack_ok=yes], [LAPACK_LIBS=])
         AC_LANG_POP([C])],
        [AC_MSG_CHECKING([for Fortran 77 LAPACK with user-supplied flags])
         AC_LANG_PUSH([Fortran 77])
         GA_F77_LAPACK_TEST()
         AC_LINK_IFELSE([], [ga_lapack_ok=yes], [LAPACK_LIBS=])
         AC_LANG_POP([Fortran 77])])
     AC_MSG_RESULT([$ga_lapack_ok])
     LIBS="$ga_save_LIBS"])

# Generic LAPACK library?
for lib in lapack lapack_rs6k; do
AS_IF([test $ga_lapack_ok = no],
    [ga_save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
     AC_CHECK_LIB([$lib], [$dgetrs],
        [ga_lapack_ok=yes; LAPACK_LIBS="-l$lib"], [], [$FLIBS])
     LIBS="$ga_save_LIBS"])
done

CPPFLAGS="$ga_save_CPPFLAGS"
LDFLAGS="$ga_save_LDFLAGS"

AC_SUBST([LAPACK_LIBS])
AC_SUBST([LAPACK_LDFLAGS])
AC_SUBST([LAPACK_CPPFLAGS])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $ga_lapack_ok = yes],
    [have_lapack=1
     $1],
    [AC_MSG_WARN([LAPACK library not found, using internal LAPACK])
     have_lapack=0
     $2])
AC_DEFINE_UNQUOTED([HAVE_LAPACK], [$have_lapack],
    [Define to 1 if using external LAPACK library])
AM_CONDITIONAL([HAVE_LAPACK], [test $ga_lapack_ok = yes])
])dnl GA_LAPACK
