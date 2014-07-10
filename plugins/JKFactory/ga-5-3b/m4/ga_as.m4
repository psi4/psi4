# GA_AS
# -----
# Certain systems may require specific assemblers (instead of $CC).
#
# Known assemblers:
# sxas  - special to NEC/NEC64
#
AC_DEFUN([GA_AS], [
AS_IF([test "x$CCAS" = x],
    [AS_CASE([$ga_cv_target], [NEC|NEC64], [CCAS=sxas])])
AS_IF([test "x$CCASFLAGS" = x],
    [AS_CASE([$ga_cv_target],
        [NEC],      [CCASFLAGS=],
        [NEC64],    [CCASFLAGS="-h size_t64"])])
])dnl
