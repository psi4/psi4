# GA_THREAD_SAFE
# --------------
# Asserts that procedures are thread safe???
# Is this specific to ARMCI???
# This was taken from older GNUmakefiles... the original doc follows
# Procedures are thread safe; should also specify the max number of threads
# via ARMCI_MAX_THREADS and thread library via THREAD_LIBRARY.  Only supported
# for SOCKETS ELAN4 OPENIB LAPI64.
AC_DEFUN([GA_THREAD_SAFE],
[AC_ARG_ENABLE([thread-safety],
    [AS_HELP_STRING([--enable-thread-safety], [**unsupported** turn on thread safety])],
    [thread_safety=yes
    AC_DEFINE([THREAD_SAFE], [1], [turn on thread safety])],
    [thread_safety=no])
AM_CONDITIONAL([THREAD_SAFE], [test x$thread_safety = xyes])
AC_ARG_VAR([THREAD_LIBRARY], [See --enable-thread-safety])
AC_SUBST([THREAD_LIBRARY])
])dnl
