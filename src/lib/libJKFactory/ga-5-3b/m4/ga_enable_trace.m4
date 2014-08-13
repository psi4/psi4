# GA_ENABLE_TRACE
# ---------------
# Whether to enable tracing. AC_DEFINEs ENABLE_TRACE.
AC_DEFUN([GA_ENABLE_TRACE],
[AC_ARG_ENABLE([trace],
    [AS_HELP_STRING([--enable-trace], [enable tracing])],
    [enable_trace=yes
    AC_DEFINE([ENABLE_TRACE], [1], [Define if tracing is enabled])],
    [enable_trace=no])
AM_CONDITIONAL([ENABLE_TRACE], [test x$enable_trace = xyes])
])dnl
