# GA_F77_UNDERSCORE
# -----------------
# If F77 requires a flag to add an underscore to external names, find it.
# Only add a single underscore regardless of whether there are already
# underscores in the name.
#
# Known flags:
#  -qextname              IBM xl
#  -qEXTNAME              IBM xl
#  -funderscoring         GNU
#  -fno-second-underscore GNU
#  -f                     absoft compiler (OSX?)
#  +ppu                   HPUX some compiler?
#
AC_DEFUN([GA_F77_UNDERSCORE],
[AC_CACHE_CHECK([for $F77 flag to add single underscore to external names],
[ga_cv_f77_underscore_flag],
[AC_LANG_PUSH([C])
AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([[void my_sub_() {}]])],
    [mv conftest.$ac_objext cfortran_test.$ac_objext
    ga_save_LIBS="$LIBS"
    LIBS="cfortran_test.$ac_objext $LIBS"
    AC_LANG_PUSH([Fortran 77])
    for ga_flag in none -qextname -qEXTNAME -funderscoring -fno-second-underscore -f +ppu ; do
        ga_save_FFLAGS="$FFLAGS"
        AS_IF([test "x$ga_flag" != xnone], [FFLAGS="$FFLAGS $ga_flag"])
        AC_LINK_IFELSE(
            [AC_LANG_CALL([], [my_sub])],
            [ga_cv_f77_underscore_flag="$ga_flag"])
        FFLAGS="$ga_save_FFLAGS"
        AS_IF([test "x$ga_cv_f77_underscore_flag" != x], [break])
    done
    AC_LANG_POP([Fortran 77])
    LIBS="$ga_save_LIBS"
    rm -f cfortran_test.$ac_objext])
AC_LANG_POP([C])])
AS_IF([test "x$ga_cv_f77_underscore_flag" != xnone],
    [AS_IF([test "x$ga_cv_f77_underscore_flag" != x],
        [FFLAGS="$FFLAGS $ga_cv_f77_underscore_flag"])])
]) # GA_F77_UNDERSCORE

# GA_F77_MAYBE_UNDERSCORING
# -------------------------
# Allow user to override the underscoring policy of their compiler.
#
# The default is to use whatever is 'natural' for the compiler, however
# GA used to be built assuming a single underscore was appended to all
# external Fortran symbols.  This allows the old behavior to be turned
# on (off by default).
#
AC_DEFUN([GA_F77_MAYBE_UNDERSCORING],
[AC_ARG_ENABLE([underscoring],
    [AS_HELP_STRING([--enable-underscoring],
        [force single underscore for all external Fortran symbols [default: use compiler defaults for external Fortran symbols]])])
AS_IF([test "x$enable_underscoring" = xyes],
    [GA_F77_UNDERSCORE])
]) # GA_F77_MAYBE_UNDERSCORING
