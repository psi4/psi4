# GA_CHECK_LIBPLOT([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------
# Check for libplot. Thus far it's only used by tcgmsg.
AC_DEFUN([GA_CHECK_LIBPLOT],
[AC_CHECK_LIB([plot], [openpl],
    [PLOTLIB=-lplot
    AC_SUBST([PLOTLIB])
    AC_DEFINE([HAVE_LIBPLOT], [1], [Defined if plot library is available])])
])dnl
