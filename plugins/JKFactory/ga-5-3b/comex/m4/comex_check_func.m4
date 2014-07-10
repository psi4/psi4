# COMEX_CHECK_FUNCS(FUNCTION..., [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------------
# Inspired by 
# http://pdh11.blogspot.com/2009/04/standard-macros-available-in-gnu.html
# but really a modified version of AC_CHECK_FUNCS.
AC_DEFUN([COMEX_CHECK_FUNCS],
[m4_map_args_w([$1], [_AH_CHECK_FUNC(], [)])]dnl
[AS_FOR([AC_func], [ac_func], [$1],
[AC_CHECK_FUNC(AC_func,
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_func), 1,
        [Define to 1 if you have AC_func, 0 if you don't]) $2],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_func), 0,
        [Define to 1 if you have AC_func, 0 if you don't]) $3], [$4])dnl])
])# COMEX_CHECK_FUNCS

# COMEX_SEARCH_LIBS(FUNCTION, LIBRARIES...,
#       [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------------
# Inspired by 
# http://pdh11.blogspot.com/2009/04/standard-macros-available-in-gnu.html
# but really a wrapped version of AC_SEARCH_LIBS.
AC_DEFUN([COMEX_SEARCH_LIBS], [
AS_VAR_PUSHDEF([HAVE_FUNC], m4_toupper(m4_translit([$1], [-.], [__])))
AC_SEARCH_LIBS([$1], [$2],
    [AC_DEFINE([HAVE_FUNC], [1], [Define to 1 if you have the `$1' function.])
     $3],
    [AC_DEFINE([HAVE_FUNC], [0], [Define to 1 if you have the `$1' function.])
     $4])
AS_VAR_POPDEF([HAVE_FUNC])
]) # COMEX_SEARCH_LIBS

