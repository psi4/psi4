# GA_CHECK_FUNCS(HEADER-FILE..., [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------------
# Inspired by 
# http://pdh11.blogspot.com/2009/04/standard-macros-available-in-gnu.html
# but really a modified version of AC_CHECK_FUNCS.
AC_DEFUN([GA_CHECK_FUNCS],
[m4_map_args_w([$1], [_AH_CHECK_FUNC(], [)])]dnl
[AS_FOR([AC_func], [ac_func], [$1],
[AC_CHECK_FUNC(AC_func,
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_func), 1,
        [Define to 1 if you have the `]AC_func[' function, 0 if you don't]) $2],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_func), 0,
        [Define to 1 if you have the `]AC_func[' function, 0 if you don't]) $3], [$4])dnl])
])# GA_CHECK_FUNCS
