# GA_PRAGMA_OPERATOR
# ------------------
# whether the C compiler understands the C99 _Pragma operator
AC_DEFUN([GA_PRAGMA_OPERATOR],
[AC_CACHE_CHECK([whether the C compiler understands the _Pragma operator],
    [ga_cv_c_pragma_operator],
    [AC_LANG_PUSH([C])
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([[_Pragma("something")]])],
        [ga_cv_c_pragma_operator=yes],
        [ga_cv_c_pragma_operator=no])
    AC_LANG_POP([C])])
AS_IF([test x$ga_cv_c_pragma_operator = xyes],
    [AC_DEFINE([HAVE_PRAGMA_OPERATOR], [1],
        [define if the C compiler understands the C99 _Pragma macro])])
AH_BOTTOM([#ifndef HAVE_PRAGMA_OPERATOR
#   define _Pragma(_pragf)
#endif])
]) #GA_PRAGMA_OPERATOR



AC_DEFUN([GA_PRAGMA_OPERATOR_OLD], [
AC_MSG_CHECKING([whether the C compiler understands the _Pragma operator])
AC_LANG_PUSH([C])
AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([[_Pragma("something")]])],
    [AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_PRAGMA_OPERATOR], [1],
        [define if the C compiler understands the C99 _Pragma macro])],
    [AC_MSG_RESULT([no])])
AC_LANG_POP([C])
]) #GA_PRAGMA_OPERATOR_OLD
