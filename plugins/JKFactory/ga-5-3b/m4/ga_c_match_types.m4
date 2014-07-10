# GA_C_POINTER_AS_INTEGER()
# -------------------------
# Size-based match between C types.
AC_DEFUN([GA_C_POINTER_AS_INTEGER],
[AC_CACHE_CHECK([for smallest C integer matching void*],
    [ga_cv_c_pointer_as_integer],
    [AS_IF(
        [test "x$ac_cv_sizeof_voidp" = "x$ac_cv_sizeof_short"],
        [ga_cv_c_pointer_as_integer=short],
        [test "x$ac_cv_sizeof_voidp" = "x$ac_cv_sizeof_int"],
        [ga_cv_c_pointer_as_integer=int],
        [test "x$ac_cv_sizeof_voidp" = "x$ac_cv_sizeof_long"],
        [ga_cv_c_pointer_as_integer=long],
        [test "x$ac_cv_sizeof_voidp" = "x$ac_cv_sizeof_long_long"],
        [ga_cv_c_pointer_as_integer="long long"],
        [AC_MSG_ERROR(
            [Could not determine smallest C integer matching void*])])])
AC_SUBST([C_POINTER_AS_INTEGER], [$ga_cv_c_pointer_as_integer])
])dnl
