# GA_ENABLE_OPT()
# ---------------
AC_DEFUN([GA_ENABLE_OPT], [
AC_ARG_ENABLE([opt],
    [AS_HELP_STRING([--disable-opt],
        [don't use hard-coded optimization flags])],
    [enable_opt=no],
    [enable_opt=yes])
])dnl
