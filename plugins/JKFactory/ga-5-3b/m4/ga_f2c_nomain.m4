# GA_F2C_NOMAIN
# ---------------
# In mixed Fortran/C code, the Fortran linker should be used. This is the
# default behavior of automake when it selects the appropriate linker.
# However, if "main" is defined in a C source file but linked by Fortran,
# there could be problems if the Fortran linker attempts to link its own main
# routine. Look for a flag to disable this behavior.
AC_DEFUN([GA_F2C_NOMAIN],
[AC_CACHE_CHECK([for flag to disable $F77 main when linking with C main],
    [ga_cv_fld_nomain],
    [AC_LANG_PUSH([C])
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[])],
        [mv conftest.$ac_objext cfortran_test.$ac_objext
        ga_save_LIBS=$LIBS
        LIBS="cfortran_test.$ac_objext $LIBS"
        AC_LANG_PUSH([Fortran 77])
        for flag in none -nofor-main -nofor_main -Mnomain -mlcmain=main; do
            ga_save_FFLAGS=$FFLAGS
            AS_IF([test "x$flag" != xnone], [FFLAGS="$FFLAGS $flag"])
            AC_LINK_IFELSE(
                [AC_LANG_SOURCE(
[[      subroutine donothing()
      end]])],
                [ga_cv_fld_nomain=$flag])
            FFLAGS=$ga_save_FFLAGS
            AS_IF([test "x$ga_cv_fld_nomain" != x], [break])
        done
        AC_LANG_POP([Fortran 77])
        LIBS=$ga_save_LIBS
        rm -f cfortran_test.$ac_objext])
    AC_LANG_POP([C])])
AS_IF([test "x$ga_cv_fld_nomain" != xnone],
    [AC_SUBST([FLD_NOMAIN], [$ga_cv_fld_nomain])])
]) # GA_F2C_NOMAIN
