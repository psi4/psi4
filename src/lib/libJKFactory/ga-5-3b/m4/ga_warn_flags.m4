# GA_WARN_FLAGS
# -------------
# Add all known warning flags to the language-specific FLAGS variable.
AC_DEFUN([GA_WARN_FLAGS], [
AS_IF([test "x$enable_warnings" = xyes], [
AC_REQUIRE([GA_COMPILER_VENDOR])
AS_VAR_PUSHDEF([vendor], [ga_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])
AS_VAR_PUSHDEF([result], [ga_cv_[]_AC_LANG_ABBREV[]_warning_flags])
AC_CACHE_CHECK([for _AC_LANG warning flags], [result], [
AC_LANG_CASE(
[C], [
AS_CASE([$vendor],
[amd],       [result="-Wall -W -Wdeclaration-after-statement"],
[borland],   [result=],
[comeau],    [result=],
[cray],      [result=],
[dec],       [result=],
[fujitsu],   [result="-Xc -pvctl,fullmsg"],
[gnu],       [result="-Wall -Wextra -Wdeclaration-after-statement -Wno-unused-parameter -pedantic -Wno-long-long -Wnested-externs -ansi"],
[hp],        [result=],
[ibm],       [result=],
[intel],     [result="-Wall"],
[kai],       [result=],
[lcc],       [result=],
[metrowerks],[result=],
[microsoft], [result=],
[pathscale], [result="-Wall -fullwarn -Wno-unused-parameter -pedantic -Wno-long-long -Wnested-externs"],
[portland],  [result="-Xc"],
[sgi],       [result=],
[sun],       [result=],
[watcom],    [result=])
],
[Fortran 77], [
result=
],
[Fortran], [
result=
],
[C++], [
result=
])
])
AC_SUBST(GA_[]_AC_LANG_PREFIX[]_WARN, [$result])
AS_VAR_POPDEF([result])
AS_VAR_POPDEF([vendor])
])
])dnl

# GA_ENABLE_WARNINGS
# ------------------
# Adds --enable-warnings.
AC_DEFUN([GA_ENABLE_WARNINGS], [
AC_ARG_ENABLE([warnings],
    [AS_HELP_STRING([--enable-warnings],
        [use compiler-specific warnings])],
    [enable_warnings=yes],
    [enable_warnings=no])
])dnl
