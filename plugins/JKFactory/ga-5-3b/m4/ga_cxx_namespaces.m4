# GA_CXX_NAMESPACES
AC_DEFUN([GA_CXX_NAMESPACES],
[AC_CACHE_CHECK([whether the compiler implements namespaces],
[ga_cv_cxx_namespaces],
[AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[namespace Outer { namespace Inner { int i = 0; }}]],
        [[using namespace Outer::Inner; return i;]])],
    [ga_cv_cxx_namespaces=yes],
    [ga_cv_cxx_namespaces=no])
AC_LANG_POP([C++])])
AS_IF([test x$ga_cv_cxx_namespaces = xyes],
    [AC_DEFINE([HAVE_NAMESPACES], [1],
        [define if the compiler implements namespaces])])
])#GA_CXX_NAMESPACES
