# GA_PROGNAME
# -----------
# Look for special global variables containing the name of the program.
AC_DEFUN([GA_PROGNAME],
[AC_CACHE_CHECK([for C global variable containing the name of the program],
    [ga_cv_progname],
    [AC_LANG_PUSH([C])
     ga_cv_progname=none
     for name in __progname program_invocation_short_name __progname_full program_invocation_name
     do
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <errno.h>
extern const char * $name;]],
[[printf("%s\n", $name);]])],
            [ga_cv_progname=$name
             break])
     done
     AC_LANG_POP([C])])
AS_IF([test x$ga_cv_progname != xnone], [val=1], [val=0])
AC_DEFINE_UNQUOTED([HAVE_PROGNAME], [$val],
    [define to 1 if the C compiler has a program name global varaible])
AC_DEFINE_UNQUOTED([PROGNAME], [$ga_cv_progname],
    [define to the name of the program name global variable])
])dnl
