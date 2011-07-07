AC_DEFUN([ACX_WITH_GOOGLE_PERF], [
  acx_with_google_perf=""
  AC_ARG_WITH([google-perf],
    [AS_HELP_STRING([--with-google-perf@<:@=Install DIR@:>@],
      [Enables use of Google fast malloc, profiler, and heap checker])],
    [
      case $withval in
      yes)
        acx_with_google_perf="yes"
      ;;
      no)
        acx_with_google_perf="no"
      ;;
      *)
        LIBS="$LIBS -L$withval/lib"
        acx_with_google_perf="$withval"
      esac
    ],
    [acx_with_google_perf="no"]
  )
  if test $acx_with_google_perf != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    if test $acx_with_libunwind != "no"; then
      AC_CHECK_LIB([tcmalloc], [malloc], [LIBS="$LIBS -ltcmalloc"], [AC_MSG_ERROR(["Unable to link with libtmalloc])])
      AC_DEFINE([MADNESS_HAS_GOOGLE_PERF], [1], [Define if using Google PerformanceTools])
    else
      AC_CHECK_LIB([tcmalloc], [malloc], [LIBS="$LIBS -ltcmalloc_minimal"], [AC_MSG_ERROR(["Unable to link with libtmalloc_minimal])])
      AC_DEFINE([MADNESS_HAS_GOOGLE_PERF_MINIMAL], [1], [Define if using Google PerformanceTools without libunwind])
    fi
    AC_LANG_RESTORE
  fi
])
