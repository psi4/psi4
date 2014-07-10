# ARMCI_STANDALONE
# ----------------
# Test whether ARMCI is configured from outside of a GA source distribution.
AC_DEFUN([ARMCI_STANDALONE],
[AC_CACHE_CHECK([whether ARMCI is configured from outside a GA distribution],
[armci_cv_standalone],
[AS_IF([test -d $ac_abs_confdir/../global],
    [armci_cv_standalone=no],
    [armci_cv_standalone=yes])])
AM_CONDITIONAL([ARMCI_STANDALONE], [test x$armci_cv_standalone = xyes])
])dnl
