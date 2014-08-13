# COMEX_DISABLE_SYS_WEAK_ALIAS
# -------------------------
# Whether to disable the test for weak aliases
AC_DEFUN([COMEX_DISABLE_SYS_WEAK_ALIAS], [
AC_ARG_ENABLE([weak],
    [AS_HELP_STRING([--disable-weak], [don't use weak symbols for profiling])],
    [],
    [enable_weak=yes])
AS_IF([test "x$comex_cv_target_base" = xCYGWIN], [enable_weak=no])
])dnl

# COMEX_SYS_WEAK_ALIAS
# -----------------
# Whether pragma weak is supported.
AC_DEFUN([COMEX_SYS_WEAK_ALIAS], [
AS_IF([test "x$enable_weak" = xyes],
    [ax_sys_weak_alias=no
     _AX_SYS_WEAK_ALIAS_PRAGMA],
    [ax_cv_sys_weak_alias_pragma=no
     AC_DEFINE([HAVE_SYS_WEAK_ALIAS_PRAGMA], [0],
        [Define this if weak aliases may be created with @%:@pragma weak])])
AM_CONDITIONAL([HAVE_SYS_WEAK_ALIAS_PRAGMA],
    [test "x$ax_cv_sys_weak_alias_pragma" = xyes])
# enable shared libs automatically if profiling using weak symbols
AS_IF([test "x$ax_cv_sys_weak_alias_pragma" = xyes],
    [AS_IF([test "x$enable_profiling" = xyes],
        [enable_shared=yes])])
])dnl
