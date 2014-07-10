# GA_CROSS_COMPILING
# ------------------
# The standard check for whether we are cross compiling is not sufficient.
# We know certain platforms are only cross compiled.  Make the fix here.
# This could be expanded later to avoid using GA_TARGET and instead perform
# a more rigorous cross compiling test case.
AC_DEFUN([GA_CROSS_COMPILING], [
AC_REQUIRE([GA_TARGET])
AC_CACHE_CHECK([whether we are cross compiling],
    [ga_cv_cross_compiling],
    [AS_IF([test "x$ga_cv_target_base" = xBGP], [cross_compiling=yes])
     ga_cv_cross_compiling=$cross_compiling])
AM_CONDITIONAL([CROSS_COMPILING], [test "x$cross_compiling" = xyes])
])dnl
