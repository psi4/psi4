# ARMCI_F77_I4
# ------------
# Force F77 to use 4-byte INTEGERs regardless of any compiler- or
# system-defaults.
AC_DEFUN([ARMCI_F77_I4], [
AC_ARG_VAR([FFLAG_INT_4],
    [Fortran 77 compiler flag to set integer size to 4 bytes])
AC_CACHE_CHECK([for INTEGER*4 size compile flag],
    [armci_cv_f77_integer_size_flag],
    [AS_IF([test x$cross_compiling = xyes],
        [_GA_F77_INTEGER_4_FLAG_CROSS([armci_cv_f77_integer_size_flag])],
        [_GA_F77_INTEGER_4_FLAG([armci_cv_f77_integer_size_flag])])])
AS_IF([test "x$armci_cv_f77_integer_size_flag" != x],
    [AS_IF([test "x$armci_cv_f77_integer_size_flag" != xnone],
        [AC_SUBST([FFLAG_INT_4], [$armci_cv_f77_integer_size_flag])])])
]) # ARMCI_F77_I4
