# _GA_F77_INTEGER_4_KNOWN_FLAGS
# -----------------------------
# These are the known flags for promoting INTEGERs to 8 bytes.
AC_DEFUN([_GA_F77_INTEGER_4_KNOWN_FLAGS],
[-fdefault-integer-4 -qintsize=4 "-integer-size 32" -CcdII4 "-s integer32" -xtypemap=integer:32 -i4 +i4])dnl

# _GA_F77_INTEGER_8_KNOWN_FLAGS
# -----------------------------
# These are the known flags for promoting INTEGERs to 8 bytes.
AC_DEFUN([_GA_F77_INTEGER_8_KNOWN_FLAGS],
[-fdefault-integer-8 -qintsize=8 "-integer-size 64" -CcdII8 "-s integer64" -xtypemap=integer:64 -i8 +i8])dnl

# _GA_F77_INTEGER_4_FLAG(VARIABLE)
# --------------------------------
# What FFLAG, if any, forces INTEGER size to be 4 bytes?
# Assign result to VARIABLE.
AC_DEFUN([_GA_F77_INTEGER_4_FLAG],
[for flag in none $FFLAG_INT _GA_F77_INTEGER_4_KNOWN_FLAGS
do
    ga_save_FFLAGS="$FFLAGS"
    AS_IF([test "x$flag" != xnone], [FFLAGS="$flag $FFLAGS"])
    sizeof_integer=0
    GA_F77_COMPUTE_SIZEOF([INTEGER], [sizeof_integer])
    FFLAGS="$ga_save_FFLAGS"
    AS_IF([test x$sizeof_integer = x4],
        [AS_TR_SH([$1])=$flag; break])
done
]) # _GA_F77_INTEGER_4_FLAG


# _GA_F77_INTEGER_8_FLAG(VARIABLE)
# --------------------------------
# What FFLAG, if any, forces INTEGER size to be 8 bytes?
# Assign result to VARIABLE.
AC_DEFUN([_GA_F77_INTEGER_8_FLAG],
[for flag in none $FFLAG_INT _GA_F77_INTEGER_8_KNOWN_FLAGS
do
    ga_save_FFLAGS="$FFLAGS"
    AS_IF([test "x$flag" != xnone], [FFLAGS="$flag $FFLAGS"])
    sizeof_integer=0
    GA_F77_COMPUTE_SIZEOF([INTEGER], [sizeof_integer])
    FFLAGS="$ga_save_FFLAGS"
    AS_IF([test x$sizeof_integer = x8],
        [AS_TR_SH([$1])=$flag; break])
done
]) # _GA_F77_INTEGER_8_FLAG


# _GA_F77_INTEGER_4_FLAG_CROSS(VARIABLE)
# --------------------------------------
# What FFLAG, if any, forces INTEGER size to be 4 bytes?
# This is safe for cross-compiling, although less accurate.
# Some compilers don't have the capability to change warnings to errors, so
# in some cases an inccorect size flag will still succeed during
# compilation. Unfortunately, there's no alternative when cross compiling.
AC_DEFUN([_GA_F77_INTEGER_4_FLAG_CROSS],
[AC_LANG_PUSH([Fortran 77])
ga_result=
ga_save_FFLAGS="$FFLAGS"
ga_save_suppress_FFLAGS="$FFLAGS"
ga_save_werror_flag=$ac_f77_werror_flag
ac_f77_werror_flag=yes
AS_IF([test "x$ga_cv_f77_suppress" != xnone],
    [ga_save_suppress_FFLAGS="$FFLAGS $ga_cv_f77_suppress"])
AS_IF([test "x$FFLAG_INT" != x],
    [FFLAGS="$ga_save_suppress_FFLAGS $FFLAG_INT"
     AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
        [ga_result=$flag])])
AS_IF([test "x$ga_result" = x],
    [for flag in _GA_F77_INTEGER_4_KNOWN_FLAGS
     do
        FFLAGS="$ga_save_suppress_FFLAGS $flag"
        AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
            [ac_ext=F
             AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
             	[ga_result=$flag; break])
             ac_ext=f])
     done])
ac_f77_werror_flag=$ga_save_werror_flag
FFLAGS="$ga_save_FFLAGS"
AS_TR_SH([$1])="$ga_result"
AC_LANG_POP([Fortran 77])
]) # _GA_F77_INTEGER_4_FLAG_CROSS


# _GA_F77_INTEGER_8_FLAG_CROSS(VARIABLE)
# --------------------------------------
# What FFLAG, if any, forces INTEGER size to be 8 bytes?
# This is safe for cross-compiling, although less accurate.
# Some compilers don't have the capability to change warnings to errors, so
# in some cases an inccorect size flag will still succeed during
# compilation. Unfortunately, there's no alternative when cross compiling.
AC_DEFUN([_GA_F77_INTEGER_8_FLAG_CROSS],
[AC_LANG_PUSH([Fortran 77])
ga_result=
ga_save_FFLAGS="$FFLAGS"
ga_save_suppress_FFLAGS="$FFLAGS"
ga_save_werror_flag=$ac_f77_werror_flag
ac_f77_werror_flag=yes
AS_IF([test "x$ga_cv_f77_suppress" != xnone],
    [ga_save_suppress_FFLAGS="$FFLAGS $ga_cv_f77_suppress"])
AS_IF([test "x$FFLAG_INT" != x],
    [FFLAGS="$ga_save_suppress_FFLAGS $FFLAG_INT"
     AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
        [ga_result=$flag])])
AS_IF([test "x$ga_result" = x],
    [for flag in _GA_F77_INTEGER_8_KNOWN_FLAGS
     do
        FFLAGS="$ga_save_suppress_FFLAGS $flag"
        AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
            [ac_ext=F
             AC_LINK_IFELSE(
[[      program main
      integer i
      end program]],
                [ga_result=$flag; break])
             ac_ext=f])
     done])
ac_f77_werror_flag=$ga_save_werror_flag
FFLAGS="$ga_save_FFLAGS"
AS_TR_SH([$1])="$ga_result"
AC_LANG_POP([Fortran 77])
]) # _GA_F77_INTEGER_8_FLAG_CROSS


# GA_F77_INTEGER_SIZE
# -------------------
# Allow the user to pick whether to use 4- or 8-byte integers.  If not
# specified, the default integer size is equivalent to sizeof(void*).
# Adds the appropriate flag to FFLAGS, if needed.
AC_DEFUN([GA_F77_INTEGER_SIZE],
[AC_ARG_VAR([FFLAG_INT],
    [Fortran 77 compiler flag to set desired integer size])
AC_ARG_ENABLE([i4],
    [AS_HELP_STRING([--enable-i4], [enable 4-byte integers [default: sizeof(void*)]])],
    [enable_i4=yes])
AC_ARG_ENABLE([i8],
    [AS_HELP_STRING([--enable-i8], [enable 8-byte integers [default: sizeof(void*)]])],
    [enable_i8=yes])
AC_LANG_PUSH([C])
AC_COMPUTE_INT([ga_f77_integer_size],
    [(long int) (sizeof (void*))],
    [AC_INCLUDES_DEFAULT()],
    [ga_f77_integer_size=0])
AC_LANG_POP([C])
AC_CACHE_CHECK([for desired Fortran INTEGER size], [ga_cv_f77_integer_size],
    [AS_IF([test x$enable_i4 = xyes],
        [AS_IF([test x$enable_i8 = xyes],
            [AC_MSG_ERROR([Cannot enable both i4 and i8])],
            [ga_cv_f77_integer_size=4])],
        [AS_IF([test x$enable_i8 = xyes],
            [ga_cv_f77_integer_size=8],
            [ga_cv_f77_integer_size=$ga_f77_integer_size])])])
# Now determine the correct compiler flag to adjust the integer size.
AC_CACHE_CHECK([for INTEGER size compile flag], [ga_cv_f77_integer_size_flag],
    [AS_CASE([$cross_compiling:$ga_cv_f77_integer_size],
        [yes:4],[_GA_F77_INTEGER_4_FLAG_CROSS([ga_cv_f77_integer_size_flag])],
        [yes:8],[_GA_F77_INTEGER_8_FLAG_CROSS([ga_cv_f77_integer_size_flag])],
        [*:4],  [_GA_F77_INTEGER_4_FLAG([ga_cv_f77_integer_size_flag])],
        [*:8],  [_GA_F77_INTEGER_8_FLAG([ga_cv_f77_integer_size_flag])])])
AS_IF([test "x$ga_cv_f77_integer_size_flag" != x],
    [AS_IF([test "x$ga_cv_f77_integer_size_flag" != xnone],
        [AC_SUBST([FFLAG_INT], [$ga_cv_f77_integer_size_flag])])])
AS_IF([test "x$ga_cv_f77_integer_size" = x8],
    [AS_IF([test "x$ga_cv_f77_integer_size_flag" = x],
        [AC_MSG_WARN([!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!])
         AC_MSG_WARN([Unable to find a flag to promote Fortran integers])
         AC_MSG_WARN([INTEGER*8 promotion is not supported])
         AC_MSG_WARN([!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!])])])
])dnl
