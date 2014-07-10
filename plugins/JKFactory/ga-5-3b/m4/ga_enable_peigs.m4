# GA_ENABLE_PEIGS
# ---------------
# Whether to enable PeIGS routines.
AC_DEFUN([GA_ENABLE_PEIGS],
[AC_ARG_ENABLE([peigs],
    [AS_HELP_STRING([--enable-peigs],
        [enable Parallel Eigensystem Solver interface])],
    [],
    [enable_peigs=no])
AS_IF([test "x$enable_peigs" = xno],
    [AC_DEFINE([ENABLE_PEIGS], [0], [Define to 1 if PeIGS is enabled])],
    [AC_DEFINE([ENABLE_PEIGS], [1], [Define to 1 if PeIGS is enabled])])
AM_CONDITIONAL([ENABLE_PEIGS], [test x$enable_peigs = xyes])
])dnl
