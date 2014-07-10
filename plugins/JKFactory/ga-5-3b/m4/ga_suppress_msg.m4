# GA_SUPPRESS_MESSAGE
# -------------------
# Determine any flags necessary to suppress informational messages from the
# compiler. Useful when setting the ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
# because any extra output (stdout and stderr) will trigger an error.
# The xlf compiler is the only case thus far.
AC_DEFUN([GA_SUPPRESS_MESSAGE],
[AC_CACHE_CHECK([for _AC_LANG flag to suppress info messages],
    [ga_cv_[]_AC_LANG_ABBREV[]_suppress],
    [ga_save_[]_AC_LANG_PREFIX[]FLAGS="$_AC_LANG_PREFIX[]FLAGS"
     ga_save_werror_flag=$ac_[]_AC_LANG_ABBREV[]_werror_flag
     ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
     for flag in none -qsuppress=cmpmsg ; do
        _AC_LANG_PREFIX[]FLAGS=$ga_save_[]_AC_LANG_PREFIX[]FLAGS
        AS_IF([test "x$flag" != xnone],
            [_AC_LANG_PREFIX[]FLAGS="$_AC_LANG_PREFIX[]FLAGS $flag"])
        AC_LINK_IFELSE([AC_LANG_PROGRAM()],
            [ga_cv_[]_AC_LANG_ABBREV[]_suppress="$flag"; break])
     done
     _AC_LANG_PREFIX[]FLAGS=$ga_save_[]_AC_LANG_PREFIX[]FLAGS
     ac_[]_AC_LANG_ABBREV[]_werror_flag=$ga_save_werror_flag
    ])
]) # GA_SUPPRESS_MESSAGE
