# GA_CHECK_PACKAGE(pkg, header, library, function, [extra-libs],
#                  [action-if-found], [action-if-not-found])
# --------------------------------------------------------------
#
AC_DEFUN([GA_CHECK_PACKAGE], [
AS_VAR_PUSHDEF([HAVE_PKG],    m4_toupper(m4_translit([HAVE_$1], [-.], [__])))
AS_VAR_PUSHDEF([PKG_LIBS],    m4_toupper(m4_translit([$1_LIBS], [-.], [__])))
AS_VAR_PUSHDEF([PKG_LDFLAGS], m4_toupper(m4_translit([$1_LDFLAGS], [-.], [__])))
AS_VAR_PUSHDEF([PKG_CPPFLAGS],m4_toupper(m4_translit([$1_CPPFLAGS], [-.], [__])))
AS_VAR_SET([PKG_LIBS],[])
AS_VAR_SET([PKG_LDFLAGS],[])
AS_VAR_SET([PKG_CPPFLAGS],[])
AC_ARG_WITH([$1],
    [AS_HELP_STRING([--with-$1[[=ARG]]],
        [specify location of $1 install and/or other flags])],
    [],
    [with_$1=yes])
AS_CASE([$with_$1],
    [yes],  [],
    [no],   [],
            [GA_ARG_PARSE(
                [with_$1],
                [PKG_LIBS],
                [PKG_LDFLAGS],
                [PKG_CPPFLAGS])])
# Check for header.
ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $PKG_CPPFLAGS"
AC_CHECK_HEADER([$2], [], [$7])
CPPFLAGS="$ga_save_CPPFLAGS"
# Check for library.
ga_save_LIBS="$LIBS"; LIBS="$PKG_LIBS $LIBS"
ga_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $PKG_LDFLAGS"
AC_SEARCH_LIBS([$4], [$3], [], [], [$5])
LIBS="$ga_save_LIBS"
LDFLAGS="$ga_save_LDFLAGS"
AS_IF([test "x$ac_cv_search_$4" != xno],
    [$6
     AC_DEFINE([HAVE_PKG], [1], [set to 1 if we have the indicated package])],
    [$7])
AS_VAR_POPDEF([HAVE_PKG])
AS_VAR_POPDEF([PKG_LIBS])
AS_VAR_POPDEF([PKG_LDFLAGS])
AS_VAR_POPDEF([PKG_CPPFLAGS])
])dnl
