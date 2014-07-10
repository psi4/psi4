# _GA_F2C_STRING(ARGLIST, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# -----------------------------------------------------------------
# Use ARGLIST in a Fortran/C mixed-language program and attempt to run it.
# The Fortran program calls a C subroutine, passing two character strings.
AC_DEFUN([_GA_F2C_STRING],
[AC_F77_FUNC([SUB])
AC_LANG_PUSH([C])
AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[#include <stdio.h>
#include <string.h>
#define LEN 80

void fix_c_string_for_f(char *s, int len);
void fix_f_string_for_c(char *s, int len);

void $SUB($1)
{
    char result[LEN];
    fix_f_string_for_c(fname, fname_len);
    fix_f_string_for_c(lname, lname_len);
    snprintf(result, LEN, "The names passed to C:  %s %s\n", fname, lname);
    fix_c_string_for_f(fname, fname_len);
    fix_c_string_for_f(lname, lname_len);
}

void fix_c_string_for_f(char *s, int len)
{
    int i;
    for (i=strlen(s); i < len; i++) {
        s[i] = ' ';
    }
}

void fix_f_string_for_c(char *s, int len)
{
    int i;
    for (i=len-1; s[i] == ' ' && i>=0; i--) {
        s[i] = '\0';
    }
}]])],
    [mv conftest.$ac_objext cfortran_test.$ac_objext
    ga_save_LIBS=$LIBS
    LIBS="cfortran_test.$ac_objext $LIBS $[]_AC_LANG_PREFIX[]LIBS"
    AC_LANG_PUSH([Fortran 77])
    AC_RUN_IFELSE([[      program main
      character(LEN=10) :: first_name
      character(LEN=15) :: last_name
      first_name = "John"
      last_name = "Doe"
      call sub(first_name, last_name)
      end]],
        [m4_default([$2], :)],
        [m4_default([$3], :)])
    AC_LANG_POP([Fortran 77])
    LIBS=$ga_save_LIBS
    rm -f cfortran_test*],
    [m4_default([$3], :)])
AC_LANG_POP([C])
rm -rf conftest*
])dnl


# GA_F2C_HIDDEN_STRING_LENGTH_CONVENTION([VALUE-IF-CROSS-COMPILING = yes])
# ---------------------------------------------------------------------------
# Determine whether a Fortran program passes the hidden string length to a
# C subroutine after all arguments versus immediately after the string.
# AC_DEFINEs F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS=<value>.
AC_DEFUN([GA_F2C_HIDDEN_STRING_LENGTH_CONVENTION],
[AC_ARG_VAR([F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS],
  [if cross compiling, set to either "yes" (default) or "no" (after string)])
AC_CACHE_CHECK([whether Fortran hidden string length convention is after args],
  [ga_cv_f2c_string_after_args],
  [AS_IF([test "x$F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS" != x],
    [ga_cv_f2c_string_after_args="$F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS"])
  AS_IF([test "x$ga_cv_f2c_string_after_args" = x],
    [AS_IF([test x$cross_compiling = xyes],
      [AS_IF([test "x$1" != x],
        [ga_cv_f2c_string_after_args="$1"],
        [ga_cv_f2c_string_after_args="yes"])])])
  AS_IF([test "x$ga_cv_f2c_string_after_args" = x],
    [_GA_F2C_STRING([char *fname, char *lname, int fname_len, int lname_len],
      [ga_cv_f2c_string_after_args="yes"])])
  AS_IF([test "x$ga_cv_f2c_string_after_args" = x],
    [_GA_F2C_STRING([char *fname, int fname_len, char *lname, int lname_len],
      [ga_cv_f2c_string_after_args="no"])])])
AS_IF([test x$cross_compiling = xyes],
  [AC_MSG_WARN([cross compiling: cannot determine f2c string convention])
   AC_MSG_WARN([default is after args (=yes) but can be changed using])
   AC_MSG_WARN([F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS=no])])
AS_IF([test "x$ga_cv_f2c_string_after_args" = "xyes"],
  [AC_DEFINE([F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS], [1],
    [whether the hidden string length comes after all other args])])
AS_IF([test "x$ga_cv_f2c_string_after_args" != x],
  [m4_default([$2], :)],
  [m4_default([$3], [AC_MSG_ERROR(
    [f2c string convention is neither after args nor after string])])])
])dnl
