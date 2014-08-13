# GA_F77_FIXED
# ------------
# Similar to AC_FC_FREEFORM, ensure that the Fortran compiler can compile 
# fixed format source code as opposed to newer free form formats.
# If necessary it will add flags to FFLAGS.
# The known flags are:
#       -ffixed-form: GNU g77, gfortran
#             -fixed: Intel compiler (ifort), Sun compiler (f95)
#            -qfixed: IBM compiler (xlf*)
#            -Mfixed: Portland Group compiler
#         -fixedform: SGI compiler
#           -f fixed: Absoft Fortran
#      +source=fixed: HP Fortran
#              -fix: Lahey/Fujitsu Fortran
#
AC_DEFUN([GA_F77_FIXED], [
AC_CACHE_CHECK([whether $F77 needs a flag to compile fixed format source],
[ga_cv_f77_fixed_flag],
[AC_LANG_PUSH([Fortran 77])
for testflag in none -ffixed-form -qfixed -Mfixed -fixedform "-f fixed" +source=fixed -fix ; do
    ga_save_FFLAGS=$FFLAGS
    AS_IF([test "x$testflag" != xnone], [FFLAGS="$FFLAGS $testflag"])
    AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[c some comment
      end program]])],
        [ga_cv_f77_fixed_flag=$testflag])
    FFLAGS=$ga_save_FFLAGS
    AS_IF([test "x$ga_cv_f77_fixed_flag" != x], [break])
done
AC_LANG_POP([Fortran 77])])
AS_IF([test "x$ga_cv_f77_fixed_flag" != xnone],
    [AS_IF([test "x$ga_cv_f77_fixed_flag" != x],
        [FFLAGS="$FFLAGS $ga_cv_f77_fixed_flag"])])
]) # GA_F77_FIXED
