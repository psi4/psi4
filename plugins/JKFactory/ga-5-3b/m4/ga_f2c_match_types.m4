# GA_F2C_MATCH_TYPES([FTYPE], [CTYPES])
# -------------------------------------
# Size-based match between Fortran and C types.
AC_DEFUN([GA_F2C_MATCH_TYPES],
[AS_VAR_PUSHDEF([ftype_size], [ga_cv_f77_sizeof_$1])
AS_VAR_PUSHDEF([ftype_cache], [ga_cv_f2c_$1])
AC_CACHE_CHECK([for C type corresponding to $1], ftype_cache,
    [m4_foreach(ctype, [$2],
        [AS_VAR_PUSHDEF([ctype_size], [ac_cv_sizeof_[]ctype])
         AS_IF([test -z "$ftype_cache"],
            [AS_IF([test "$ctype_size" = "$ftype_size"],
                [ftype_cache="ctype"])])
         AS_VAR_POPDEF([ctype_size])])])
AS_IF([test "x$ftype_cache" = x],
    [AC_MSG_ERROR([Could not determine C type matching Fortran $1])])
AC_SUBST(AS_TR_SH([F2C_$1_C_TYPE]), [$ftype_cache])
AS_VAR_POPDEF([ftype_cache])
AS_VAR_POPDEF([ftype_size])
])dnl
