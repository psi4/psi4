# GA_FUNCTION
# -----------
# Define FUNCTION_NAME to either __func__ or __FUNCTION__ appropriately.
# If all else fails, #define FUNCTION_NAME <blank>.
AC_DEFUN([GA_FUNCTION],
[AC_CACHE_CHECK([for preprocessor symbol for function name],
[ga_cv_cpp_function],
[AS_IF([test x$ga_cv_cpp_function = x],
    [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[extern int printf(const char *format, ...);]],
            [[printf("__func__ = %s\n", __func__);]])],
        [ga_cv_cpp_function=__func__])])
AS_IF([test x$ga_cv_cpp_function = x],
    [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[extern int printf(const char *format, ...);]],
            [[printf("__FUNCTION__ = %s\n", __FUNCTION__);]])],
        [ga_cv_cpp_function=__FUNCTION__])])
AS_IF([test x$ga_cv_cpp_function = x],
    [ga_cv_cpp_function='"Unknown"'])])
AC_DEFINE_UNQUOTED([FUNCTION_NAME], [$ga_cv_cpp_function],
    [CPP symbol for function name, if available])
])# GA_FUNCTION
