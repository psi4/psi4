# GA_F77_MISMATCH_TYPE
# --------------------
# Determine F77 compiler flags, if any, that are required for compiling code
# which uses inconsistent data types for certain function calls.  Nearly all
# of our Fortran routines are written in C and some are not type safe.
# The known flags are:
#   -mismatch_all: nagfor
#   -dusty: nagfor
#
AC_DEFUN([GA_F77_MISMATCH_TYPE], [
AC_CACHE_CHECK([whether $F77 needs a flag to compile inconsistent types],
[ga_cv_f77_mismatch_type_flag],
[AC_LANG_PUSH([Fortran 77])
for testflag in none "-mismatch_all -dusty"
do
    ga_save_FFLAGS="$FFLAGS"
    AS_IF([test "x$testflag" != xnone], [FFLAGS="$FFLAGS $testflag"])
    AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE(
[[      integer function the_test ()
      implicit none
      logical foo
      external foo
      character*1       byte_mb
      integer           int_mb
      the_test = 0
      if (foo(byte_mb)) return
      if (foo( int_mb)) return
      the_test = 1
      return
      end]])],
        [ga_cv_f77_mismatch_type_flag="$testflag"])
    FFLAGS="$ga_save_FFLAGS"
    AS_IF([test "x$ga_cv_f77_mismatch_type_flag" != x], [break])
done
AC_LANG_POP([Fortran 77])])
AS_IF([test "x$ga_cv_f77_mismatch_type_flag" != xnone],
    [AS_IF([test "x$ga_cv_f77_mismatch_type_flag" != x],
        [FFLAGS="$FFLAGS $ga_cv_f77_mismatch_type_flag"])])
]) # GA_F77_MISMATCH_TYPE
