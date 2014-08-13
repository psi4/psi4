# COMEX_CHECK_HEADERS(HEADER-FILE...,
#                  [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND], [INCLUDES])
# ----------------------------------------------------------------------
# Inspired by 
# http://pdh11.blogspot.com/2009/04/standard-macros-available-in-gnu.html
# but really a modified version of AC_CHECK_HEADERS.
AC_DEFUN([COMEX_CHECK_HEADERS],
[m4_map_args_w([$1], [_AH_CHECK_HEADER(], [)])]dnl
[AS_FOR([AC_header], [ac_header], [$1],
[AC_CHECK_HEADER(AC_header,
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_header), 1,
        [Define to 1 if you have AC_header, 0 if you don't]) $2],
    [AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_]AC_header), 0,
        [Define to 1 if you have AC_header, 0 if you don't]) $3], [$4])dnl])
])# COMEX_CHECK_HEADERS
