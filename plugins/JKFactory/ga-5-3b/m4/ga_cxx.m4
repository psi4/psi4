# GA_CXX_ENABLE()
# ---------------
# Whether to enable C++ bindings.
AC_DEFUN([GA_CXX_ENABLE],
[AC_ARG_ENABLE([cxx],
    [AS_HELP_STRING([--enable-cxx],[build C++ interface])],
    [enable_cxx=yes],
    [enable_cxx=no])
AM_CONDITIONAL([CXX_BINDINGS],[test x$enable_cxx = xyes])
])# GA_CXX_ENABLE
