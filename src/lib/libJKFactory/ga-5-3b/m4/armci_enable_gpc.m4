# ARMCI_ENABLE_GPC
# ----------------
# Whether to enable GPC calls in ARMCI.
AC_DEFUN([ARMCI_ENABLE_GPC], [
AC_ARG_ENABLE([gpc],
    [AS_HELP_STRING([--enable-gpc], [enable GPC calls])],
    [enable_gpc=yes],
    [enable_gpc=no])
AS_IF([test $enable_gpc = yes],
    [AC_DEFINE([ARMCI_ENABLE_GPC_CALLS], [1],
        [Define to 1 if GPC calls are enabled])],
    [AC_DEFINE([ARMCI_ENABLE_GPC_CALLS], [0],
        [Define to 1 if GPC calls are enabled])])
AM_CONDITIONAL([ARMCI_ENABLE_GPC_CALLS], [test x$enable_gpc = xyes])
])dnl
