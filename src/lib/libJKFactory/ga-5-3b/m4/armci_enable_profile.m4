# ARMCI_ENABLE_PROFILE
# --------------------
# Whether to enable profiling. AC_DEFINEs ARMCI_PROFILE.
AC_DEFUN([ARMCI_ENABLE_PROFILING], [
AC_ARG_ENABLE([profiling],
    [AS_HELP_STRING([--enable-profiling], [enable profiling])],
    [],
    [enable_profiling=no])
AS_IF([test "x$enable_profiling" = xyes],
    [AC_DEFINE([ARMCI_PROFILE], [1], [Define if profiling is enabled])])
AM_CONDITIONAL([ENABLE_PROFILING], [test "x$enable_profiling" = xyes])
])dnl
