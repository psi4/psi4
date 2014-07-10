# COMEX_AR
# -----
# Libtool doesn't advertise AR nor AR_FLAGS in case the user wishes to
# override them. Further, certain systems require a different archiver.
# RANLIB may also be affected.
# Use this prior to LT_INIT.
#
# Known archivers:
# ar    - all known systems
# sxar  - special to NEC/NEC64
#
AC_DEFUN([COMEX_AR], [
AC_ARG_VAR([AR], [archiver used by libtool (default: ar)])
AC_ARG_VAR([AR_FLAGS], [archiver flags used by libtool (default: cru)])
AC_ARG_VAR([RANLIB], [generates index to archive (default: ranlib)])
AS_IF([test "x$AR" = x],
    [AS_CASE([$comex_cv_target], [NEC|NEC64], [AR=sxar])])
AS_IF([test "x$RANLIB" = x],
    [AS_CASE([$comex_cv_target], [NEC|NEC64], [RANLIB=true])])
])dnl
