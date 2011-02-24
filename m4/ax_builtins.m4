# SYNOPSIS
#
#   AX_BUILTINS()
#
# DESCRIPTION
#
#   This macro checks compiler support for builtin functions.
#   If the checked functions are supported a macro is defined.

AC_DEFUN([AX_BUILTINS], [
    AC_PREREQ(2.59)

    # see if the C++ compiler supports __builtin_expect
    AC_LANG_PUSH(C++)
    AC_CACHE_CHECK([if $CXX supports __builtin_expect],
        [_cv_cxx_supports___builtin_expect],
        [AC_TRY_LINK([],
          [void *ptr = (void*) 0;
           if (__builtin_expect (ptr != (void*) 0, 1)) return 0;],
          [_cv_cxx_supports___builtin_expect="yes"],
          [_cv_cxx_supports___builtin_expect="no"])])
    if test "$_cv_cxx_supports___builtin_expect" = "yes" ; then
        have_builtin_expect=1
    else
        have_builtin_expect=0
    fi
    AC_DEFINE_UNQUOTED([HAVE_BUILTIN_EXPECT], [$have_builtin_expect],
          [Whether C++ compiler supports __builtin_expect])
    AC_LANG_POP(C++)

    # see if the C compiler supports __builtin_prefetch
    AC_LANG_PUSH(C++)
    AC_CACHE_CHECK([if $CXX supports __builtin_prefetch],
        [_cv_cxx_supports___builtin_prefetch],
        [AC_TRY_LINK([],
          [int ptr;
           __builtin_prefetch(&ptr,0,0);],
          [_cv_cxx_supports___builtin_prefetch="yes"],
          [_cv_cxx_supports___builtin_prefetch="no"])])
    if test "$_cv_cxx_supports___builtin_prefetch" = "yes" ; then
        have_builtin_prefetch=1
    else
        have_builtin_prefetch=0
    fi
    AC_DEFINE_UNQUOTED([HAVE_BUILTIN_PREFETCH], [$have_builtin_prefetch],
          [Whether C++ compiler supports __builtin_prefetch])
    AC_LANG_POP(C++)

]) dnl AX_BUILTINS
