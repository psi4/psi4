# MA_LONG_DOUBLE_TYPEDEF
# ----------------------
# AC_SUBST a proper long double if one does not exist.
AC_DEFUN([MA_LONG_DOUBLE_TYPEDEF],
[AC_REQUIRE([AC_TYPE_LONG_DOUBLE])
AS_IF([test $ac_cv_type_long_double = yes],
    [ma_long_double="long double"],
    [ma_long_double="struct {double dummy[2];}"])
AC_SUBST([MA_LONG_DOUBLE], [$ma_long_double])
]) # MA_LONG_DOUBLE_TYPEDEF
