# GA_F77_FLUSH
# ------------
# Check whether the function "flush" is available. If not, try various flags
# in order to link to it successfully.  Otherwise look for the "forflush"
# function.  Add any appropriate flags to FFLAGS as needed and define
# F77_FLUSH to the name of the appropriate function.
# If flush routine is available set HAVE_F77_FLUSH to 1, otherwise 0.
AC_DEFUN([GA_F77_FLUSH], [
AC_CACHE_CHECK([for $F77 flush routine],
[ga_cv_f77_flush],
[AC_LANG_PUSH([Fortran 77])
for testflag in none -qextname=flush ; do
    ga_save_FFLAGS=$FFLAGS
    AS_IF([test "x$testflag" != xnone], [FFLAGS="$FFLAGS $testflag"])
    AC_LINK_IFELSE([AC_LANG_SOURCE(
[[      call flush(13)
      end program]])],
        [ga_cv_f77_flush_flag=$testflag])
    FFLAGS=$ga_save_FFLAGS
    AS_IF([test "x$ga_cv_f77_flush_flag" != x], [break])
done
AS_IF([test "x$ga_cv_f77_flush_flag" != x], [ga_cv_f77_flush=flush],
    [AC_LINK_IFELSE([AC_LANG_SOURCE(
[[      call forflush(13)
      end program]])],
        [ga_cv_f77_flush=forflush])])
AC_LANG_POP([Fortran 77])])
AS_IF([test "x$ga_cv_f77_flush" != x],
    [AC_DEFINE_UNQUOTED([F77_FLUSH], [$ga_cv_f77_flush],
        [Name of F77 flush routine])
     AC_DEFINE([HAVE_F77_FLUSH], [1],
        [whether F77 flush routine is available])],
    [AC_MSG_WARN([Could not determine name of $F77 flush routine])
     AC_DEFINE([HAVE_F77_FLUSH], [0],
        [whether F77 flush routine is available])])
AS_IF([test "x$ga_cv_f77_flush" = xflush],
    [AS_IF([test "x$ga_cv_f77_flush_flag" != xnone],
        [FFLAGS="$FFLAGS $ga_cv_f77_flush_flag"])])
]) # GA_F77_FLUSH
