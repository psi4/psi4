# GA_ENABLE_PROFILING
# -------------------
# Whether to enable profiling. AC_DEFINEs ENABLE_PROFILING.
AC_DEFUN([GA_ENABLE_PROFILING],
[AC_ARG_ENABLE([profiling],
    [AS_HELP_STRING([--enable-profiling], [enable profiling])],
    [],
    [enable_profiling=no])
AM_CONDITIONAL([ENABLE_PROFILING], [test x$enable_profiling = xyes])
AS_IF([test "x$enable_profiling" = xyes],
    [AC_DEFINE([ENABLE_PROFILING], [1], [set to 1 if profiling is enabled])],
    [AC_DEFINE([ENABLE_PROFILING], [0], [set to 1 if profiling is enabled])])
])dnl
