###############################################################################
# Look for Global Arrays
###############################################################################
AC_DEFUN([GA_REQUIRE], [
AC_CACHE_CHECK([whether this package is part of GA distribution],
    [ga_cv_libga_la_found],
    [AS_IF([test -f "../libga.la"],
        [ga_cv_libga_la_found=yes],
        [ga_cv_libga_la_found=no])])
AS_IF([test "x$ga_cv_libga_la_found" = xyes],
    [GA_LIBS="../libga.la"
     GA_CPPFLAGS="-I$srcdir/../global/src -I../global/src -I$srcdir/../ma -I../ma -I../gaf2c -I../armci/gaf2c -I$srcdir/../armci/src/include"
     GA_FLIBS=`echo '@FLIBS@' | ../config.status --file=-`],
    [ga_save_PATH="$PATH"
     AS_IF([test -d $with_ga], [PATH="$with_ga:$PATH"])
     AS_IF([test -d $with_ga/bin], [PATH="$with_ga/bin:$PATH"])
     AC_PATH_PROG([GA_CONFIG], [ga-config])
     PATH="$ga_save_PATH"
     AS_IF([test "x$GA_CONFIG" != x],
        [GA_CPPFLAGS=`$GA_CONFIG --cppflags`
         GA_LDFLAGS=`$GA_CONFIG --ldflags`
         GA_LIBS=`$GA_CONFIG --libs`
         GA_FLIBS=`$GA_CONFIG --flibs`],
        [GA_CHECK_PACKAGE([ga], [ga.h], [ga], [GA_Initialize], [$FLIBS],
            [GA_LIBS="-lga"],
            [AS_UNSET([ac_cv_search_GA_Initialize])
             GA_CHECK_PACKAGE([ga], [ga.h], [global], [GA_Initialize],
                [-lma -larmci -llinalg $FLIBS],
                [GA_LIBS="-lglobal -lma -larmci -llinalg"],
                [AC_MSG_FAILURE([Could not locate Global Arrays 5.x nor 4.x])])])])])
AC_SUBST([GA_CPPFLAGS])
AC_SUBST([GA_LDFLAGS])
AC_SUBST([GA_LIBS])
AC_SUBST([GA_FLIBS])
])dnl
