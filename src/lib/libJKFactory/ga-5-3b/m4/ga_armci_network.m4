# _GA_ARMCI_NETWORK_WITH(KEY, DESCRIPTION)
# --------------------------------------------------
# A helper macro for generating all of the AC_ARG_WITHs.
# Also may establish value of ga_armci_network.
# Counts how many armci networks were specified by user.
AC_DEFUN([_GA_ARMCI_NETWORK_WITH], [
AC_ARG_WITH([$1],
    [AS_HELP_STRING([--with-$1[[=ARG]]], [select armci network as $2])])
AS_VAR_PUSHDEF([KEY],      m4_toupper(m4_translit([$1],      [-.], [__])))
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
dnl Can't have AM_CONDITIONAL here in case configure must find armci network
dnl without user intervention.
dnl AM_CONDITIONAL([ARMCI_NETWORK_]KEY, [test "x$with_key" != x])
AS_IF([test "x$with_key" != x],
    [GA_ARG_PARSE([with_key], [ARMCI_NETWORK_LIBS], [ARMCI_NETWORK_LDFLAGS],
                  [ARMCI_NETWORK_CPPFLAGS])])
AS_IF([test "x$with_key" != xno && test "x$with_key" != x],
    [ga_armci_network=KEY
     AS_VAR_ARITH([armci_network_count], [$armci_network_count + 1])])
AS_VAR_POPDEF([KEY])
AS_VAR_POPDEF([with_key])
])dnl

# _GA_ARMCI_NETWORK_WARN(KEY)
# ---------------------------
# Helper macro for idicating value of armci network arguments.
AC_DEFUN([_GA_ARMCI_NETWORK_WARN], [
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
AS_IF([test "x$with_key" != x && test "x$with_key" != xno],
    [AC_MSG_WARN([--with-$1=$with_key])])
AS_VAR_POPDEF([with_key])
])dnl

# _GA_ARMCI_NETWORK_AM_CONDITIONAL(KEY)
#--------------------------------------
# Helper macro for generating all AM_CONDITIONALs.
AC_DEFUN([_GA_ARMCI_NETWORK_AM_CONDITIONAL], [
AS_VAR_PUSHDEF([KEY],      m4_toupper(m4_translit([$1],      [-.], [__])))
AS_VAR_PUSHDEF([with_key],            m4_translit([with_$1], [-.], [__]))
AM_CONDITIONAL([ARMCI_NETWORK_]KEY,
    [test "x$with_key" != x && test "x$with_key" != xno])
AS_VAR_POPDEF([KEY])
AS_VAR_POPDEF([with_key])
])dnl

# _GA_ARMCI_NETWORK_ARMCI([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_ARMCI], [
AC_MSG_NOTICE([searching for external ARMCI...])
happy=yes
CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
LIBS="$LIBS $GA_MP_LIBS"
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([armci.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([ARMCI_Init], [armci], [], [happy=no])
     AS_CASE([$ac_cv_search_ARMCI_Init],
            ["none required"], [],
            [no], [],
            [# add missing lib to ARMCI_NETWORK_LIBS if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*$ac_cv_search_ARMCI_Init*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS $ac_cv_search_ARMCI_Init"])])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([armci_group_comm], [armci])
     AS_IF([test "x$ac_cv_search_armci_group_comm" != xno],
        [ac_cv_search_armci_group_comm=1],
        [ac_cv_search_armci_group_comm=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_GROUP_COMM],
        [$ac_cv_search_armci_group_comm],
        [set to 1 if ARMCI has armci_group_comm function])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_MEMBER([ARMCI_Group.comm], [], [], [[#include <armci.h>]])
     AS_IF([test "x$ac_cv_member_ARMCI_Group_comm" != xno],
        [ac_cv_member_ARMCI_Group_comm=1],
        [ac_cv_member_ARMCI_Group_comm=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_GROUP_COMM_MEMBER],
        [$ac_cv_member_ARMCI_Group_comm],
        [set to 1 if ARMCI has ARMCI_Group.comm member])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([ARMCI_Initialized], [armci])
     AS_IF([test "x$ac_cv_search_ARMCI_Initialized" != xno],
        [ac_cv_search_ARMCI_Initialized=1],
        [ac_cv_search_ARMCI_Initialized=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_INITIALIZED],
        [$ac_cv_search_ARMCI_Initialized],
        [set to 1 if ARMCI has ARMCI_Initialized function])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([armci_stride_info_init], [armci])
     AS_IF([test "x$ac_cv_search_armci_stride_info_init" != xno],
        [ac_cv_search_armci_stride_info_init=1],
        [ac_cv_search_armci_stride_info_init=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_STRIDE_INFO_INIT],
        [$ac_cv_search_armci_stride_info_init],
        [set to 1 if ARMCI has armci_stride_info_init function])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([armci_notify], [armci])
     AS_IF([test "x$ac_cv_search_armci_notify" != xno],
        [ac_cv_search_armci_notify=1],
        [ac_cv_search_armci_notify=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_NOTIFY],
        [$ac_cv_search_armci_notify],
        [set to 1 if ARMCI has armci_notify function])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([armci_msg_init], [armci])
     AS_IF([test "x$ac_cv_search_armci_msg_init" != xno],
        [ac_cv_search_armci_msg_init=1],
        [ac_cv_search_armci_msg_init=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_MSG_INIT],
        [$ac_cv_search_armci_msg_init],
        [set to 1 if ARMCI has armci_msg_init function])
    ])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([armci_msg_finalize], [armci])
     AS_IF([test "x$ac_cv_search_armci_msg_finalize" != xno],
        [ac_cv_search_armci_msg_finalize=1],
        [ac_cv_search_armci_msg_finalize=0])
     AC_DEFINE_UNQUOTED([HAVE_ARMCI_MSG_FINALIZE],
        [$ac_cv_search_armci_msg_finalize],
        [set to 1 if ARMCI has armci_msg_finalize function])
    ])
AM_CONDITIONAL([HAVE_ARMCI_GROUP_COMM],
   [test "x$ac_cv_search_armci_group_comm" = x1])
AM_CONDITIONAL([HAVE_ARMCI_GROUP_COMM_MEMBER],
   [test "x$ac_cv_member_ARMCI_Group_comm" = x1])
AM_CONDITIONAL([HAVE_ARMCI_INITIALIZED],
   [test "x$ac_cv_search_ARMCI_Initialized" = x1])
AM_CONDITIONAL([HAVE_ARMCI_STRIDE_INFO_INIT],
   [test "x$ac_cv_search_armci_stride_info_init" = x1])
AM_CONDITIONAL([HAVE_ARMCI_NOTIFY],
   [test "x$ac_cv_search_armci_notify" = x1])
AM_CONDITIONAL([HAVE_ARMCI_MSG_INIT],
   [test "x$ac_cv_search_armci_msg_init" = x1])
AM_CONDITIONAL([HAVE_ARMCI_MSG_FINALIZE],
   [test "x$ac_cv_search_armci_msg_finalize" = x1])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=ARMCI; with_armci=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_BGML([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_BGML], [
AC_MSG_NOTICE([searching for BGML...])
happy=yes
dnl AS_IF([test "x$happy" = xyes],
dnl     [AS_IF([test -d /bgl/BlueLight/ppcfloor/bglsys], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([BGLML_memcpy], [msglayer.rts], [], [happy=no],
        [-lrts.rts -ldevices.rts])
     AS_CASE([$ac_cv_search_BGLML_memcpy],
        ["none required"], [],
        [no], [],
        [# add msglayer.rts to ARMCI_NETWORK_LIBS if not there
         AS_CASE([$ARMCI_NETWORK_LIBS],
                 [*msglayer.rts*], [],
                 [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -lmsglayer.rts"])
         # add extra lib rts.rts to ARMCI_NETWORK_LIBS if not there
         AS_CASE([$ARMCI_NETWORK_LIBS],
                 [*rts.rts*], [],
                 [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -lrts.rts"])
         # add extra lib devices.rts to ARMCI_NETWORK_LIBS if not there
         AS_CASE([$ARMCI_NETWORK_LIBS],
                 [*devices.rts*], [],
                 [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -ldevices.rts"])])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=BGML; with_bgml=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_CRAY_SHMEM([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_CRAY_SHMEM], [
AC_MSG_NOTICE([searching for CRAY_SHMEM...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([mpp/shmem.h], [],
        [AC_CHECK_HEADER([shmem.h], [], [happy=no])])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([shmem_init], [sma], [], [happy=no])
     AS_CASE([$ac_cv_search_shmem_init],
        ["none required"], [],
        [no], [],
        [# add sma to ARMCI_NETWORK_LIBS if not there
         AS_CASE([$ARMCI_NETWORK_LIBS],
                 [*sma*], [],
                 [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -lsma"])])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=CRAY_SHMEM; with_cray_shmem=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_DCMF([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_DCMF], [
AC_MSG_NOTICE([searching for DCMF...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([dcmf.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([DCMF_Messager_initialize], [dcmf.cnk],
        [], [happy=no], [-ldcmfcoll.cnk -lSPI.cna -lrt])
     AS_CASE([$ac_cv_search_DCMF_Messager_initialize],
            ["none required"], [],
            [no], [],
            [# add dcmf.cnk to ARMCI_NETWORK_LIBS if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*dcmf.cnk*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -ldcmf.cnk"])
             # add extra lib dcmfcoll.cnk if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*dcmfcoll.cnk*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -ldcmfcoll.cnk"])
             # add extra lib SPI.cna if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*SPI.cna*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -lSPI.cna"])
             # add extra lib rt if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*rt*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS -lrt"])])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=DCMF; with_dcmf=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_LAPI([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ----------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_LAPI], [
AC_MSG_NOTICE([searching for LAPI...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([lapi.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([LAPI_Init], [lapi_r lapi], [], [happy=no])
     AS_CASE([$ac_cv_search_LAPI_Init],
            ["none required"], [],
            [no], [],
            [# add missing lib to ARMCI_NETWORK_LIBS if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*$ac_cv_search_LAPI_Init*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS $ac_cv_search_LAPI_Init"])])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=LAPI; with_lapi=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_MPI_TS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_MPI_TS], [
AC_MSG_NOTICE([searching for MPI_TS...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=MPI_TS; with_mpi_ts=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_MPI_MT([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_MPI_MT], [
AC_MSG_NOTICE([searching for MPI_MT...])
happy=yes
CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
LIBS="$LIBS $GA_MP_LIBS"
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([mpi.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([MPI_Init_thread], [mpi mpich.cnk mpich.rts],
        [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=MPI_MT; with_mpi_mt=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_MPI_SPAWN([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_MPI_SPAWN], [
AC_MSG_NOTICE([searching for MPI_SPAWN...])
happy=yes
CPPFLAGS="$CPPFLAGS $GA_MP_CPPFLAGS"
LDFLAGS="$LDFLAGS $GA_MP_LDFLAGS"
LIBS="$LIBS $GA_MP_LIBS"
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([mpi.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([MPI_Comm_spawn_multiple], [mpi mpich.cnk mpich.rts],
        [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=MPI_SPAWN; with_mpi_spawn=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_OFA([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_OFA], [
AC_MSG_NOTICE([searching for OFA...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([infiniband/verbs.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([ibv_open_device], [ibverbs], [], [happy=no])
     AS_CASE([$ac_cv_search_ibv_open_device],
        ["none required"], [],
        [no], [],
        [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS $ac_cv_search_ibv_open_device"])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=OFA; with_ofa=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_OPENIB([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_OPENIB], [
AC_MSG_NOTICE([searching for OPENIB...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([infiniband/verbs.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([ibv_open_device], [ibverbs], [], [happy=no])
     AS_CASE([$ac_cv_search_ibv_open_device],
        ["none required"], [],
        [no], [],
        [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS $ac_cv_search_ibv_open_device"])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=OPENIB; with_openib=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_PORTALS4([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_PORTALS4], [
AC_MSG_NOTICE([searching for PORTALS4...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=PORTALS4; with_portals4=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_PORTALS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -------------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_PORTALS], [
AC_MSG_NOTICE([searching for PORTALS...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([portals/portals3.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([portals/nal.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([PtlInit], [portals], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=PORTALS; with_portals=yes; $1],
    [$2])
])dnl

# _GA_ARMCI_NETWORK_DMAPP([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------
AC_DEFUN([_GA_ARMCI_NETWORK_DMAPP], [
AC_MSG_NOTICE([searching for DMAPP...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=DMAPP; with_dmapp=yes; $1],
    [$2])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([gethugepagesize], [hugetlbfs])
     AS_CASE([$ac_cv_search_gethugepagesize],
            ["none required"], [],
            [no], [],
            [# add missing lib to ARMCI_NETWORK_LIBS if not there
             AS_CASE([$ARMCI_NETWORK_LIBS],
                     [*$ac_cv_search_gethugepagesize*], [],
                     [ARMCI_NETWORK_LIBS="$ARMCI_NETWORK_LIBS $ac_cv_search_gethugepagesize"])])])
])dnl

# _GA_ARMCI_NETWORK_GEMINI([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
# TODO when gemini headers and libraries become available, fix this
AC_DEFUN([_GA_ARMCI_NETWORK_GEMINI], [
AC_MSG_NOTICE([searching for GEMINI...])
happy=yes
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([numatoolkit.h], [], [happy=no], [
AC_INCLUDES_DEFAULT
#include <mpi.h>])])
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([NTK_Init], [numatoolkit], [], [happy=no])])
# CPPFLAGS must have CRAY_UGNI before looking for the next headers.
gemini_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS -DCRAY_UGNI"
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([onesided.h], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [AC_CHECK_HEADER([gni.h], [], [happy=no])])
CPPFLAGS="$gemini_save_CPPFLAGS"
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([gniInit], [onesided], [], [happy=no])])
AS_IF([test "x$happy" = xyes],
    [ga_armci_network=GEMINI; with_gemini=yes; $1],
    [$2])
# check for a function introduced in libonesided/1.5
# we purposefully abuse the ac_cv_search_onesided_mem_htflush value
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([onesided_mem_htflush], [onesided])
     AS_IF([test "x$ac_cv_search_onesided_mem_htflush" != xno],
        [ac_cv_search_onesided_mem_htflush=1],
        [ac_cv_search_onesided_mem_htflush=0])
     AC_DEFINE_UNQUOTED([HAVE_ONESIDED_MEM_HTFLUSH],
        [$ac_cv_search_onesided_mem_htflush],
        [set to 1 if libonesided has onesided_mem_htflush (added in v1.5)])
    ])
# check for a function introduced in libonesided/1.6
# we purposefully abuse the ac_cv_search_onesided_fadd value
AS_IF([test "x$happy" = xyes],
    [AC_SEARCH_LIBS([onesided_fadd], [onesided])
     AS_IF([test "x$ac_cv_search_onesided_fadd" != xno],
        [ac_cv_search_onesided_fadd=1],
        [ac_cv_search_onesided_fadd=0])
     AC_DEFINE_UNQUOTED([HAVE_ONESIDED_FADD],
        [$ac_cv_search_onesided_fadd],
        [set to 1 if libonesided has onesided_fadd (added in v1.6)])
    ])
])dnl

# GA_ARMCI_NETWORK
# ----------------
# This macro allows user to choose the armci network but also allows the
# network to be tested for automatically.
AC_DEFUN([GA_ARMCI_NETWORK], [
# Clear the variables we will be using, just in case.
ga_armci_network=
ARMCI_NETWORK_LIBS=
ARMCI_NETWORK_LDFLAGS=
ARMCI_NETWORK_CPPFLAGS=
AC_ARG_ENABLE([autodetect],
    [AS_HELP_STRING([--enable-autodetect],
        [attempt to locate ARMCI_NETWORK besides SOCKETS])])
# First, all of the "--with" stuff is taken care of.
armci_network_count=0
_GA_ARMCI_NETWORK_WITH([armci],     [external; path to external ARMCI library])
_GA_ARMCI_NETWORK_WITH([bgml],      [IBM BG/L])
_GA_ARMCI_NETWORK_WITH([cray-shmem],[Cray XT shmem])
_GA_ARMCI_NETWORK_WITH([dcmf],      [IBM BG/P Deep Computing Message Framework])
_GA_ARMCI_NETWORK_WITH([dmapp],     [(Comex) Cray DMAPP])
_GA_ARMCI_NETWORK_WITH([gemini],    [Cray XE Gemini using libonesided])
_GA_ARMCI_NETWORK_WITH([lapi],      [IBM LAPI])
_GA_ARMCI_NETWORK_WITH([mpi-mt],    [MPI-2 multi-threading])
_GA_ARMCI_NETWORK_WITH([mpi-spawn], [MPI-2 dynamic process mgmt])
_GA_ARMCI_NETWORK_WITH([mpi-ts],    [(Comex) MPI-1 two-sided])
_GA_ARMCI_NETWORK_WITH([ofa],       [(Comex) Infiniband OpenIB])
_GA_ARMCI_NETWORK_WITH([openib],    [Infiniband OpenIB])
_GA_ARMCI_NETWORK_WITH([portals4],  [(Comex) Portals4])
_GA_ARMCI_NETWORK_WITH([portals],   [Cray XT portals])
_GA_ARMCI_NETWORK_WITH([sockets],   [Ethernet TCP/IP])
# Temporarily add ARMCI_NETWORK_CPPFLAGS to CPPFLAGS.
ga_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $ARMCI_NETWORK_CPPFLAGS"
# Temporarily add ARMCI_NETWORK_LDFLAGS to LDFLAGS.
ga_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $ARMCI_NETWORK_LDFLAGS"
# Temporarily add ARMCI_NETWORK_LIBS to LIBS.
ga_save_LIBS="$LIBS"; LIBS="$ARMCI_NETWORK_LIBS $LIBS"
AS_IF([test "x$enable_autodetect" = xyes],
    [AC_MSG_NOTICE([searching for ARMCI_NETWORK...])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_bgml" != xno],
        [_GA_ARMCI_NETWORK_BGML()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_cray_shmem" != xno],
        [_GA_ARMCI_NETWORK_CRAY_SHMEM()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_dcmf" != xno],
        [_GA_ARMCI_NETWORK_DCMF()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_lapi" != xno],
        [_GA_ARMCI_NETWORK_LAPI()])
dnl     AS_IF([test "x$ga_armci_network" = x && test "x$with_mpi_ts" != xno],
dnl         [_GA_ARMCI_NETWORK_MPI_TS()])
dnl     AS_IF([test "x$ga_armci_network" = x && test "x$with_mpi_mt" != xno],
dnl         [_GA_ARMCI_NETWORK_MPI_MT()])
dnl     AS_IF([test "x$ga_armci_network" = x && test "x$with_mpi_spawn" != xno],
dnl         [_GA_ARMCI_NETWORK_MPI_SPAWN()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_ofa" != xno],
        [_GA_ARMCI_NETWORK_OFA()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_openib" != xno],
        [_GA_ARMCI_NETWORK_OPENIB()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_portals4" != xno],
        [_GA_ARMCI_NETWORK_PORTALS4()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_portals" != xno],
        [_GA_ARMCI_NETWORK_PORTALS()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_dmapp" != xno],
        [_GA_ARMCI_NETWORK_DMAPP()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_gemini" != xno],
        [_GA_ARMCI_NETWORK_GEMINI()])
     AS_IF([test "x$ga_armci_network" = x && test "x$with_armci" != xno],
        [_GA_ARMCI_NETWORK_ARMCI()])
     AS_IF([test "x$ga_armci_network" = x],
        [AC_MSG_WARN([!!!])
         AC_MSG_WARN([No ARMCI_NETWORK detected, defaulting to SOCKETS])
         AC_MSG_WARN([!!!])
         ga_armci_network=SOCKETS; with_sockets=yes])],
    [# Not autodetecting
     # Check whether multiple armci networks were selected by user.
     AS_CASE([$armci_network_count],
        [0], [AC_MSG_WARN([No ARMCI_NETWORK specified, defaulting to SOCKETS])
              ga_armci_network=SOCKETS; with_sockets=yes],
        [1], [AS_IF([test "x$ga_armci_network" = xARMCI],
                 [_GA_ARMCI_NETWORK_ARMCI([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=ARMCI failed])])])
              AS_IF([test "x$ga_armci_network" = xBGML],
                 [_GA_ARMCI_NETWORK_BGML([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=BGML failed])])])
              AS_IF([test "x$ga_armci_network" = xCRAY_SHMEM],
                 [_GA_ARMCI_NETWORK_CRAY_SHMEM([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=CRAY_SHMEM failed])])])
              AS_IF([test "x$ga_armci_network" = xDCMF],
                 [_GA_ARMCI_NETWORK_DCMF([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=DCMF failed])])])
              AS_IF([test "x$ga_armci_network" = xDMAPP],
                 [_GA_ARMCI_NETWORK_DMAPP([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=DMAPP failed])])])
              AS_IF([test "x$ga_armci_network" = xLAPI],
                 [_GA_ARMCI_NETWORK_LAPI([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=LAPI failed])])])
              AS_IF([test "x$ga_armci_network" = xMPI_TS],
                 [_GA_ARMCI_NETWORK_MPI_TS([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=MPI_TS failed])])])
              AS_IF([test "x$ga_armci_network" = xMPI_MT],
                 [_GA_ARMCI_NETWORK_MPI_MT([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=MPI_MT failed])])])
              AS_IF([test "x$ga_armci_network" = xMPI_SPAWN],
                 [_GA_ARMCI_NETWORK_MPI_SPAWN([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=MPI_SPAWN failed])])])
              AS_IF([test "x$ga_armci_network" = xOFA],
                 [_GA_ARMCI_NETWORK_OFA([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=OFA failed])])])
              AS_IF([test "x$ga_armci_network" = xOPENIB],
                 [_GA_ARMCI_NETWORK_OPENIB([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=OPENIB failed])])])
              AS_IF([test "x$ga_armci_network" = xPORTALS4],
                 [_GA_ARMCI_NETWORK_PORTALS4([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=PORTALS4 failed])])])
              AS_IF([test "x$ga_armci_network" = xPORTALS],
                 [_GA_ARMCI_NETWORK_PORTALS([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=PORTALS failed])])])
              AS_IF([test "x$ga_armci_network" = xGEMINI],
                 [_GA_ARMCI_NETWORK_GEMINI([],
                    [AC_MSG_ERROR([test for ARMCI_NETWORK=GEMINI failed])])])
             ],
        [AC_MSG_WARN([too many armci networks specified: $armci_network_count])
         AC_MSG_WARN([the following were specified:])
         _GA_ARMCI_NETWORK_WARN([armci])
         _GA_ARMCI_NETWORK_WARN([bgml])
         _GA_ARMCI_NETWORK_WARN([cray-shmem])
         _GA_ARMCI_NETWORK_WARN([dcmf])
         _GA_ARMCI_NETWORK_WARN([dmapp])
         _GA_ARMCI_NETWORK_WARN([lapi])
         _GA_ARMCI_NETWORK_WARN([mpi-ts])
         _GA_ARMCI_NETWORK_WARN([mpi-mt])
         _GA_ARMCI_NETWORK_WARN([mpi-spawn])
         _GA_ARMCI_NETWORK_WARN([ofa])
         _GA_ARMCI_NETWORK_WARN([openib])
         _GA_ARMCI_NETWORK_WARN([portals4])
         _GA_ARMCI_NETWORK_WARN([portals])
         _GA_ARMCI_NETWORK_WARN([gemini])
         _GA_ARMCI_NETWORK_WARN([sockets])
         AC_MSG_ERROR([please select only one armci network])])])
# Remove ARMCI_NETWORK_CPPFLAGS from CPPFLAGS.
CPPFLAGS="$ga_save_CPPFLAGS"
# Remove ARMCI_NETWORK_LDFLAGS from LDFLAGS.
LDFLAGS="$ga_save_LDFLAGS"
# Remove ARMCI_NETWORK_LIBS from LIBS.
LIBS="$ga_save_LIBS"
_GA_ARMCI_NETWORK_AM_CONDITIONAL([armci])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([bgml])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([cray-shmem])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([dcmf])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([dmapp])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([lapi])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([mpi-ts])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([mpi-mt])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([mpi-spawn])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([ofa])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([openib])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([gemini])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([portals4])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([portals])
_GA_ARMCI_NETWORK_AM_CONDITIONAL([sockets])
AC_SUBST([ARMCI_NETWORK_LDFLAGS])
AC_SUBST([ARMCI_NETWORK_LIBS])
AC_SUBST([ARMCI_NETWORK_CPPFLAGS])

# TODO
AM_CONDITIONAL([DCMF_VER_2],   [test x != x])  # always false
AM_CONDITIONAL([DCMF_VER_0_2], [test x != x]) # always false
AM_CONDITIONAL([DCMF_VER_0_3], [test x = x]) # always true

# permanent hack
AS_CASE([$ga_armci_network],
[DMAPP],    [ARMCI_SRC_DIR=comex],
[GEMINI],   [ARMCI_SRC_DIR=src-gemini],
[MPI_MT],   [ARMCI_SRC_DIR=src],
[MPI_TS],   [ARMCI_SRC_DIR=comex],
[OFA],      [ARMCI_SRC_DIR=comex],
[OPENIB],   [ARMCI_SRC_DIR=src],
[PORTALS4], [ARMCI_SRC_DIR=comex],
[PORTALS],  [ARMCI_SRC_DIR=src-portals],
[GEMINI],   [ARMCI_SRC_DIR=src-gemini],
            [ARMCI_SRC_DIR=src])
AC_SUBST([ARMCI_SRC_DIR])
AM_CONDITIONAL([ARMCI_SRC_DIR_PORTALS], [test "x$ARMCI_SRC_DIR" = "xsrc-portals"])
AM_CONDITIONAL([ARMCI_SRC_DIR_GEMINI],  [test "x$ARMCI_SRC_DIR" = "xsrc-gemini"])
AM_CONDITIONAL([ARMCI_SRC_DIR_COMEX],   [test "x$ARMCI_SRC_DIR" = "xcomex"])
AM_CONDITIONAL([ARMCI_SRC_DIR_SRC],     [test "x$ARMCI_SRC_DIR" = "xsrc"])

# tcgmsg5 requires this
AS_IF([test x$ga_armci_network = xLAPI],
[AC_DEFINE([NOTIFY_SENDER], [1],
    [this was defined unconditionally when using LAPI for tcgmsg 5])
AC_DEFINE([LAPI], [1], [tcgmsg 5 requires this when using LAPI])
])

ga_cray_xt_networks=no
AS_IF([test x$ga_armci_network = xPORTALS], [ga_cray_xt_networks=yes])
AS_IF([test x$ga_armci_network = xCRAY_SHMEM], [ga_cray_xt_networks=yes])
AM_CONDITIONAL([CRAY_XT_NETWORKS], [test x$ga_cray_xt_networks = xyes])

ga_cv_sysv_hack=no
# Only perform this hack for ARMCI build.
AS_IF([test "x$ARMCI_TOP_BUILDDIR" != x], [
    AS_IF([test x$ga_cv_sysv = xno],
        [AS_CASE([$ga_armci_network],
            [BGML|DCMF|PORTALS|GEMINI], [ga_cv_sysv_hack=no],
                [ga_cv_sysv_hack=yes])],
        [ga_cv_sysv_hack=yes])
AS_IF([test x$ga_cv_sysv_hack = xyes],
    [AC_DEFINE([SYSV], [1],
        [Defined if we want this system to use SYSV shared memory])])
])
AM_CONDITIONAL([SYSV], [test x$ga_cv_sysv_hack = xyes])

# if not using external armci library, the following functions are always available
AS_IF([test "x$ga_armci_network" != xARMCI],
    [AC_DEFINE([HAVE_ARMCI_GROUP_COMM], [1], [])
     AC_DEFINE([HAVE_ARMCI_INITIALIZED], [1], [])
     AC_DEFINE([HAVE_ARMCI_NOTIFY], [1], [])
     AC_DEFINE([HAVE_ARMCI_MSG_INIT], [1], [])
     AC_DEFINE([HAVE_ARMCI_MSG_FINALIZE], [1], [])])
AM_CONDITIONAL([HAVE_ARMCI_GROUP_COMM_MEMBER],
   [test "x$ac_cv_member_ARMCI_Group_comm" = x1])
AM_CONDITIONAL([HAVE_ARMCI_GROUP_COMM],  [test "x$ga_armci_network" != xARMCI])
AM_CONDITIONAL([HAVE_ARMCI_INITIALIZED], [test "x$ga_armci_network" != xARMCI])
AM_CONDITIONAL([HAVE_ARMCI_NOTIFY],      [test "x$ga_armci_network" != xARMCI])
AM_CONDITIONAL([HAVE_ARMCI_MSG_INIT],    [test "x$ga_armci_network" != xARMCI])
AM_CONDITIONAL([HAVE_ARMCI_MSG_FINALIZE],[test "x$ga_armci_network" != xARMCI])
# the armci iterators only available in the conglomerate sources
AS_CASE([$ga_armci_network],
    [ARMCI|GEMINI|PORTALS], [],
    [AC_DEFINE([HAVE_ARMCI_STRIDE_INFO_INIT], [1], [])])
AM_CONDITIONAL([HAVE_ARMCI_STRIDE_INFO_INIT],
    [test "x$ga_armci_network" != xARMCI && test "x$ga_armci_network" != xGEMINI && test "x$ga_armci_network" != xPORTALS])

# ugly hack for working around NWChem memory requirements
# and MPI startup verus the 'classic' ARMCI startup
delay_tcgmsg_mpi_startup=1
AS_CASE([$ga_armci_network],
[ARMCI],        [delay_tcgmsg_mpi_startup=0],
[BGML],         [delay_tcgmsg_mpi_startup=1],
[CRAY_SHMEM],   [delay_tcgmsg_mpi_startup=1],
[DCMF],         [delay_tcgmsg_mpi_startup=1],
[DMAPP],        [delay_tcgmsg_mpi_startup=0],
[LAPI],         [delay_tcgmsg_mpi_startup=1],
[MPI_TS],       [delay_tcgmsg_mpi_startup=0],
[MPI_SPAWN],    [delay_tcgmsg_mpi_startup=1],
[OFA],          [delay_tcgmsg_mpi_startup=0],
[OPENIB],       [delay_tcgmsg_mpi_startup=1],
[GEMINI],       [delay_tcgmsg_mpi_startup=1],
[PORTALS4],     [delay_tcgmsg_mpi_startup=0],
[PORTALS],      [delay_tcgmsg_mpi_startup=1],
[SOCKETS],      [delay_tcgmsg_mpi_startup=1])
AC_DEFINE_UNQUOTED([NEED_DELAY_TCGMSG_MPI_STARTUP],
    [$delay_tcgmsg_mpi_startup],
    [whether to wait until the last moment to call ARMCI_Init() in TCGMSG-MPI])
])dnl
