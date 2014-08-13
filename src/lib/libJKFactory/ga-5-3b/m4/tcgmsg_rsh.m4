# TCGMSG_REMOTE_SHELL
# -------------------
# tcgmsg requires a remote shell.
AC_DEFUN([TCGMSG_REMOTE_SHELL],
[AC_MSG_NOTICE([checking for remote shell])
AC_PATH_PROGS([ga_cv_path_rsh], [rsh remsh ssh], [not found])
AS_IF([test "x$ga_cv_path_rsh" = "xnot found"],
    [AC_MSG_ERROR([Could not find remote shell for use by TCGMSG])])
AC_DEFINE_UNQUOTED([TCGMSG_RSH], ["$ga_cv_path_rsh"], [remote shell for TCGMSG])
])dnl
