# _COMEX_NETWORK_WITH(KEY, DESCRIPTION)
# --------------------------------------------------
# A helper macro for generating all of the AC_ARG_WITHs.
# Also may establish value of comex_network.
# Counts how many comex networks were specified by user.
AC_DEFUN([_COMEX_NETWORK_WITH], [
AC_ARG_WITH([$1],
    [AS_HELP_STRING([--with-$1[[=ARG]]], [select comex network as $2])])
AS_VAR_PUSHDEF([KEY],      m4_toupper(m4_translit([$1],      [-.], [__])))
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
dnl Can't have AM_CONDITIONAL here in case configure must find comex network
dnl without user intervention.
dnl AM_CONDITIONAL([COMEX_NETWORK_]KEY, [test "x$with_key" != x])
AS_IF([test "x$with_key" != x],
    [COMEX_ARG_PARSE([with_key], [COMEX_NETWORK_LIBS], [COMEX_NETWORK_LDFLAGS],
                  [COMEX_NETWORK_CPPFLAGS])])
AS_IF([test "x$with_key" != xno && test "x$with_key" != x],
    [comex_network=KEY
     AS_VAR_ARITH([comex_network_count], [$comex_network_count + 1])])
AS_VAR_POPDEF([KEY])
AS_VAR_POPDEF([with_key])
])dnl

# _COMEX_NETWORK_WARN(KEY)
# ---------------------------
# Helper macro for idicating value of comex network arguments.
AC_DEFUN([_COMEX_NETWORK_WARN], [
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
AS_IF([test "x$with_key" != x && test "x$with_key" != xno],
    [AC_MSG_WARN([--with-$1=$with_key])])
AS_VAR_POPDEF([with_key])
])dnl

# _COMEX_NETWORK_AC_DEFINE(KEY)
#--------------------------------------
# Helper macro for generating all AC_DEFINEs.
AC_DEFUN([_COMEX_NETWORK_AC_DEFINE], [
AS_VAR_PUSHDEF([KEY],      m4_toupper(m4_translit([$1],      [-.], [__])))
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
AS_IF([test "x$with_key" != x && test "x$with_key" != xno],
    [AC_DEFINE([COMEX_NETWORK_]KEY, [1], [Define to 1 if the network is ]KEY)],
    [AC_DEFINE([COMEX_NETWORK_]KEY, [0], [Define to 1 if the network is ]KEY)])
AS_VAR_POPDEF([KEY])
AS_VAR_POPDEF([with_key])
])dnl

# _COMEX_NETWORK_AM_CONDITIONAL(KEY)
#--------------------------------------
# Helper macro for generating all AM_CONDITIONALs.
AC_DEFUN([_COMEX_NETWORK_AM_CONDITIONAL], [
AS_VAR_PUSHDEF([KEY],      m4_toupper(m4_translit([$1],      [-.], [__])))
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
AM_CONDITIONAL([COMEX_NETWORK_]KEY,
    [test "x$with_key" != x && test "x$with_key" != xno])
AS_VAR_POPDEF([KEY])
AS_VAR_POPDEF([with_key])
])dnl

# _COMEX_NETWORK_MPI_TS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([_COMEX_NETWORK_MPI_TS], [
AC_MSG_NOTICE([searching for MPI_TS...])
happy=yes
CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
LDFLAGS="$LDFLAGS $MPI_LDFLAGS"
LIBS="$LIBS $MPI_LIBS"
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([mpi.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([MPI_Init], [], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [comex_network=MPI_TS; with_mpi_ts=yes; $1],
    [$2])
])dnl

# _COMEX_NETWORK_OFA([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
AC_DEFUN([_COMEX_NETWORK_OFA], [
AC_MSG_NOTICE([searching for OFA...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([infiniband/verbs.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([ibv_open_device], [ibverbs], [], [happy=no])
     AS_CASE([$ac_cv_search_ibv_open_device],
        ["none required"], [],
        [no], [],
        [COMEX_NETWORK_LIBS="$COMEX_NETWORK_LIBS $ac_cv_search_ibv_open_device"])])
AS_IF([test "x$happy" = xyes],
    [comex_network=OFA; with_ofa=yes; $1],
    [$2])
])dnl

# _COMEX_NETWORK_PORTALS4([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -------------------------------------------------------------------
AC_DEFUN([_COMEX_NETWORK_PORTALS4], [
AC_MSG_NOTICE([searching for PORTALS4...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([portals4.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([PtlInit], [portals4], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([PtlFini], [portals4], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [comex_network=PORTALS4; with_portals4=yes; $1],
    [$2])
])dnl

# _COMEX_NETWORK_DMAPP([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
# TODO when dmapp headers and libraries become available, fix this
AC_DEFUN([_COMEX_NETWORK_DMAPP], [
AC_MSG_NOTICE([searching for DMAPP...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([dmapp.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([gethugepagesize], [hugetlbfs],
        [AC_DEFINE([HAVE_LIBHUGETLBFS], [1],
               [Define to 1 if you have the `hugetlbfs' library.])],
        [AC_DEFINE([HAVE_LIBHUGETLBFS], [0],
               [Define to 1 if you have the `hugetlbfs' library.])])
     AS_CASE([$ac_cv_search_gethugepagesize],
            ["none required"], [],
            [no], [],
            [# add missing lib to COMEX_NETWORK_LIBS if not there
             AS_CASE([$COMEX_NETWORK_LIBS],
                     [*$ac_cv_search_gethugepagesize*], [],
                     [COMEX_NETWORK_LIBS="$COMEX_NETWORK_LIBS $ac_cv_search_gethugepagesize"])])
     AC_CHECK_TYPES([dmapp_lock_desc_t], [], [], [[#include <dmapp.h>]])
     AC_CHECK_TYPES([dmapp_lock_handle_t], [], [], [[#include <dmapp.h>]])
    ])
AS_IF([test "x$happy" = xyes],
    [comex_network=DMAPP; with_dmapp=yes; $1],
    [$2])
])dnl

# COMEX_NETWORK_SETUP
# -------------------
# This macro allows user to choose the comex network but also allows the
# network to be tested for automatically.
AC_DEFUN([COMEX_NETWORK_SETUP], [
# Clear the variables we will be using, just in case.
comex_network=
COMEX_NETWORK_LIBS=
COMEX_NETWORK_LDFLAGS=
COMEX_NETWORK_CPPFLAGS=
AC_ARG_ENABLE([autodetect],
    [AS_HELP_STRING([--enable-autodetect],
        [attempt to locate COMEX_NETWORK besides MPI two-sided])])
# First, all of the "--with" stuff is taken care of.
comex_network_count=0
_COMEX_NETWORK_WITH([mpi-ts],    [MPI-1 two-sided])
_COMEX_NETWORK_WITH([ofa],       [Infiniband OpenIB])
_COMEX_NETWORK_WITH([portals4],  [Portals4])
_COMEX_NETWORK_WITH([dmapp],     [Cray DMAPP])
# Temporarily add COMEX_NETWORK_CPPFLAGS to CPPFLAGS.
comex_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $COMEX_NETWORK_CPPFLAGS"
# Temporarily add COMEX_NETWORK_LDFLAGS to LDFLAGS.
comex_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $COMEX_NETWORK_LDFLAGS"
# Temporarily add COMEX_NETWORK_LIBS to LIBS.
comex_save_LIBS="$LIBS"; LIBS="$COMEX_NETWORK_LIBS $LIBS"
AS_IF([test "x$enable_autodetect" = xyes],
    [AC_MSG_NOTICE([searching for COMEX_NETWORK...])
     AS_IF([test "x$comex_network" = x && test "x$with_ofa" != xno],
        [_COMEX_NETWORK_OFA()])
     AS_IF([test "x$comex_network" = x && test "x$with_portals4" != xno],
        [_COMEX_NETWORK_PORTALS4()])
     AS_IF([test "x$comex_network" = x && test "x$with_dmapp" != xno],
        [_COMEX_NETWORK_DMAPP()])
     AS_IF([test "x$comex_network" = x],
        [AC_MSG_WARN([!!!])
         AC_MSG_WARN([No COMEX_NETWORK detected, defaulting to MPI_TS])
         AC_MSG_WARN([!!!])
         comex_network=MPI_TS; with_mpi_ts=yes])],
    [# Not autodetecting
     # Check whether multiple comex networks were selected by user.
     AS_CASE([$comex_network_count],
        [0], [AC_MSG_WARN([No COMEX_NETWORK specified, defaulting to MPI_TS])
              comex_network=MPI_TS; with_mpi_ts=yes],
        [1], [AS_IF([test "x$comex_network" = xMPI_TS],
                 [_COMEX_NETWORK_MPI_TS([],
                    [AC_MSG_ERROR([test for COMEX_NETWORK=MPI_TS failed])])])
              AS_IF([test "x$comex_network" = xOFA],
                 [_COMEX_NETWORK_OFA([],
                    [AC_MSG_ERROR([test for COMEX_NETWORK=OFA failed])])])
              AS_IF([test "x$comex_network" = xPORTALS4],
                 [_COMEX_NETWORK_PORTALS4([],
                    [AC_MSG_ERROR([test for COMEX_NETWORK=PORTALS4 failed])])])
              AS_IF([test "x$comex_network" = xDMAPP],
                 [_COMEX_NETWORK_DMAPP([],
                    [AC_MSG_ERROR([test for COMEX_NETWORK=DMAPP failed])])])
             ],
        [AC_MSG_WARN([too many comex networks specified: $comex_network_count])
         AC_MSG_WARN([the following were specified:])
         _COMEX_NETWORK_WARN([mpi-ts])
         _COMEX_NETWORK_WARN([ofa])
         _COMEX_NETWORK_WARN([portals4])
         _COMEX_NETWORK_WARN([dmapp])
         AC_MSG_ERROR([please select only one comex network])])])
# Remove COMEX_NETWORK_CPPFLAGS from CPPFLAGS.
CPPFLAGS="$comex_save_CPPFLAGS"
# Remove COMEX_NETWORK_LDFLAGS from LDFLAGS.
LDFLAGS="$comex_save_LDFLAGS"
# Remove COMEX_NETWORK_LIBS from LIBS.
LIBS="$comex_save_LIBS"
_COMEX_NETWORK_AM_CONDITIONAL([mpi-ts])
_COMEX_NETWORK_AM_CONDITIONAL([ofa])
_COMEX_NETWORK_AM_CONDITIONAL([portals4])
_COMEX_NETWORK_AM_CONDITIONAL([dmapp])
_COMEX_NETWORK_AC_DEFINE([mpi-ts])
_COMEX_NETWORK_AC_DEFINE([ofa])
_COMEX_NETWORK_AC_DEFINE([portals4])
_COMEX_NETWORK_AC_DEFINE([dmapp])
AC_SUBST([COMEX_NETWORK_LDFLAGS])
AC_SUBST([COMEX_NETWORK_LIBS])
AC_SUBST([COMEX_NETWORK_CPPFLAGS])
])dnl
