# use AX_BOOST_BASE first, then try to find boost headers using other means
AC_DEFUN([GA_BOOST],
[AX_BOOST_BASE([1.35.0], [boost_ok=yes], [boost_ok=no])
 # could not find boost install
 # try $ac_boost_path or $BOOST_ROOT as CPPFLAG
 AS_IF([test "x$boost_ok" = xno],
    [AC_CACHE_CHECK(
        [for non-installed and non-staged boost headers],
        [ga_cv_boost_headers],
        [gfutex_save_CPPFLAGS="$CPPFLAGS"
         CPPFLAGS="$CPPFLAGS -I$ac_boost_path"
         AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
            [[@%:@include <boost/version.hpp>]],
            [[#if BOOST_VERSION >= $WANT_BOOST_VERSION
            // Everything is okay
            #else
            #  error Boost version is too old
            #endif]])],
            [ga_cv_boost_headers="$ac_boost_path"])
         CPPFLAGS="$gfutex_save_CPPFLAGS"
         AS_IF([test "x$ga_cv_boost_headers" = x],
            [gfutex_save_CPPFLAGS="$CPPFLAGS"
             CPPFLAGS="$CPPFLAGS -I$BOOST_ROOT"
             AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                [[@%:@include <boost/version.hpp>]],
                [[#if BOOST_VERSION >= $WANT_BOOST_VERSION
                // Everything is okay
                #else
                #  error Boost version is too old
                #endif]])],
                [ga_cv_boost_headers="$BOOST_ROOT"])
             CPPFLAGS="$gfutex_save_CPPFLAGS"])
         AS_IF([test "x$ga_cv_boost_headers" = x],
            [ga_cv_boost_headers="no"])])
     AS_IF([test "x$ga_cv_boost_headers" = xno],
        [AC_MSG_ERROR([could not locate boost headers])],
        [BOOST_CPPFLAGS="-I$ga_cv_boost_headers"])
     AC_SUBST([BOOST_CPPFLAGS])])
])
