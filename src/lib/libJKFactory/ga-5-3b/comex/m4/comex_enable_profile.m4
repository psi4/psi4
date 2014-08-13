# COMEX_ENABLE_PROFILING_ARMCI
# ----------------------------
# Whether to enable profiling for ARMCI.
# AC_DEFINEs COMEX_PROFILING_ARMCI.
AC_DEFUN([COMEX_ENABLE_PROFILING_ARMCI], [
AC_ARG_ENABLE([profiling-armci],
    [AS_HELP_STRING([--enable-profiling-armci],
                    [enable profiling for ARMCI])],
    [],
    [enable_profiling_armci=no])
AS_IF([test "x$enable_profiling_armci" = xyes],
    [AC_DEFINE([COMEX_PROFILING_ARMCI], [1],
               [Define if ARMCI profiling is enabled])])
AM_CONDITIONAL([ENABLE_PROFILING_ARMCI],
               [test "x$enable_profiling_armci" = xyes])
])dnl

# COMEX_ENABLE_PROFILING_COMEX
# ----------------------------
# Whether to enable profiling for COMEX.
# AC_DEFINEs COMEX_PROFILING_COMEX.
AC_DEFUN([COMEX_ENABLE_PROFILING_COMEX], [
AC_ARG_ENABLE([profiling-comex],
    [AS_HELP_STRING([--enable-profiling-comex],
                    [enable profiling for ComEx])],
    [],
    [enable_profiling_comex=no])
AS_IF([test "x$enable_profiling_comex" = xyes],
    [AC_DEFINE([COMEX_PROFILING_COMEX], [1],
               [Define if ComEx profiling is enabled])])
AM_CONDITIONAL([ENABLE_PROFILING_COMEX],
               [test "x$enable_profiling_comex" = xyes])
])dnl
