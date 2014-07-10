# GA_WITH_HELP
# ------------
# Using undocumented features, add some text to our --with-PACKAGE help to
# avoid repetition.
AC_DEFUN([GA_WITH_HELP],
[AC_ARG_WITH([PACKAGE], [AS_HELP_STRING([--with-PACKAGE[[=ARG]]],
    [for most of the external software packages, ARG can be one or more whitespace-separated directories, linker or preprocessor directives; for example, --with-PACKAGE="/path/to/PACKAGE -lmylib -I/mydir"])])])
