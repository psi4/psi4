# ARMCI_ENABLE_NB_NONCONT
# -----------------------
# Not sure what this is for.
AC_DEFUN([ARMCI_ENABLE_NB_NONCONT],
[AC_ARG_ENABLE([nb_noncont],
    [AS_HELP_STRING([--enable-nb-noncont], [TODO])],
    [enable_nb_noncont=yes
    AC_DEFINE([NB_NONCONT], [1], [TODO])],
    [enable_nb_noncont=no])
AM_CONDITIONAL([NB_NONCONT], [test x$enable_nb_noncont = xyes])
])dnl
