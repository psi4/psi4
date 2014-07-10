# GA_F77_DISABLE()
# ----------------
# Whether to disable all Fortran code.
AC_DEFUN([GA_F77_DISABLE],
[AC_ARG_ENABLE([f77],
    [AS_HELP_STRING([--disable-f77], [disable Fortran code])],
    [],
    [enable_f77=yes])
])# GA_F77_DISABLE

# GA_F77_DISABLE_RESULTS()
# ------------------------
# This used to be part of the the GA_F77_DISABLE macro, but it turns out the
# fortran compiler may be bogus and therefore we should disable it after
# GA_F77_DISABLE has been run. These AC_DEFINEs and AM_CONDITIONALs should be
# set at the last possible moment when the final value for enable_f77 has been
# set.
AC_DEFUN([GA_F77_DISABLE_RESULTS], [
AS_IF([test "x$enable_f77" = xyes],
    [AC_DEFINE([NOFORT],     [0], [Define to 1 if not using Fortran])
     AC_DEFINE([ENABLE_F77], [1], [Define to 1 if using Fortran])],
    [AC_DEFINE([NOFORT],     [1], [Define to 1 if not using Fortran])
     AC_DEFINE([ENABLE_F77], [0], [Define to 1 if using Fortran])])
AM_CONDITIONAL([NOFORT],     [test "x$enable_f77" = xno])
AM_CONDITIONAL([ENABLE_F77], [test "x$enable_f77" = xyes])
])# GA_F77_DISABLE_RESULTS
