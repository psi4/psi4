# ARMCI_ENABLE_GROUP
# ------------------
# Not sure what this is for.
AC_DEFUN([ARMCI_ENABLE_GROUP],
[AC_ARG_ENABLE([armci_group],
    [AS_HELP_STRING([--enable-armci-group], [TODO])],
    [enable_armci_group=yes
    AC_DEFINE([ARMCI_GROUP], [1], [TODO])],
    [enable_armci_group=no])
AM_CONDITIONAL([ARMCI_GROUP], [test x$enable_armci_group = xyes])
])dnl
