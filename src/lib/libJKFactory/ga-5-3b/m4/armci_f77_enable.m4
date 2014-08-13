# ARMCI_F77_ENABLE()
# ------------------
# Whether to enable Fortran code within ARMCI.
AC_DEFUN([ARMCI_F77_ENABLE],
[AC_ARG_ENABLE([f77],
    [AS_HELP_STRING([--enable-f77], [enable Fortran code])],
    [],
    [enable_f77=no])
AS_IF([test "x$enable_f77" = xyes],
    [AC_DEFINE([NOFORT],     [0], [Define to 1 if not using Fortran])
     AC_DEFINE([ENABLE_F77], [1], [Define to 1 if using Fortran])],
    [AC_DEFINE([NOFORT],     [1], [Define to 1 if not using Fortran])
     AC_DEFINE([ENABLE_F77], [0], [Define to 1 if using Fortran])])
AM_CONDITIONAL([NOFORT],     [test "x$enable_f77" = xno])
AM_CONDITIONAL([ENABLE_F77], [test "x$enable_f77" = xyes])
])# ARMCI_F77_ENABLE
