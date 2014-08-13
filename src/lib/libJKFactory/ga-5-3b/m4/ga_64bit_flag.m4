# GA_64BIT_FLAG
# -------------
# Figure out whether the compiler needs a special flag to generate 64bit code.
# We verify by testing the size of void*. We can't use AC_CHECK_SIZEOF
# directly since it uses the cache. Instead we copied code from autoconf...
#
# Known flags:
#  -m64     GNU
#  -q64     IBM
#  +DD64    HPUX
#  +DA2.0W  HPUX (obsolete form of +DD64)
#  -64      SGI TFP, not sure, others might be -mips64, -align64?
#
AC_DEFUN([GA_64BIT_FLAG],
[AC_CACHE_CHECK([for flag to indicate 64-bits], [ga_cv_64bit_flag],
[AC_LANG_PUSH([C])
for flag in none -m64 -q64 +DD64 +DA2.0W -64 ; do
    ga_save_CFLAGS=$CFLAGS
    AS_IF([test "x$flag" != xnone], [CFLAGS="$CFLAGS $flag"])
    ga_sizeof_voidp=0
    AC_COMPUTE_INT([ga_sizeof_voidp],
        [(long int) (sizeof (void*))],
        [AC_INCLUDES_DEFAULT()],
        [ga_sizeof_voidp=0])
    CFLAGS=$ga_save_CFLAGS
    AS_IF([test x$ga_sizeof_voidp = x8],
        [ga_cv_64bit_flag=$flag; break])
done
AC_LANG_POP([C])])
AS_IF([test x$ga_cv_64bit_flag != xnone],
    [CFLAGS="$CFLAGS $ga_cv_64bit_flag"
     FFLAGS="$FFLAGS $ga_cv_64bit_flag"])
]) # GA_64BIT_FLAG
