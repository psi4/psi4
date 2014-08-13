# GA_TARGET()
# -----------
# Attempt to determine the old TARGET variable automatically.
# Deprecated TARGETs:
# CRAY-SV1
# cray-sv2
# CRAY-T3E
# CRAY-YMP
# CYGNUS
# DECOSF
# HITACHI
# INTERIX
# SGI
# SGI64
# SGI_N32
# SGITFP
AC_DEFUN([GA_TARGET],
[# AH_TEMPLATE for all known TARGETs
AH_TEMPLATE([BGL],          [Define to 1 on BlueGene/L systems])
AH_TEMPLATE([BGP],          [Define to 1 on BlueGene/P systems])
AH_TEMPLATE([CATAMOUNT],    [Define to 1 on Cray XT systems using Catamount])
AH_TEMPLATE([CRAY_SV1],     [Define to 1 on Cray SV1 systems])
AH_TEMPLATE([CRAY_SV2],     [Define to 1 on Cray SV2 systems])
AH_TEMPLATE([CRAY_T3E],     [Define to 1 on Cray T3E systems])
AH_TEMPLATE([CRAY_XT],      [Define to 1 on Cray XT systems])
AH_TEMPLATE([CRAY_YMP],     [Define to 1 on Cray YMP systems])
AH_TEMPLATE([CYGNUS],       [Define to 1 on Cygnus systems])
AH_TEMPLATE([CYGWIN],       [Define to 1 on Cygwin systems])
AH_TEMPLATE([DECOSF],       [Define to 1 on DEC OSF])
AH_TEMPLATE([FUJITSU_VPP],  [Define to 1 on fujitsu systems])
AH_TEMPLATE([FUJITSU_VPP64],[Define to 1 on fujitsu systems])
AH_TEMPLATE([HITACHI],      [Define to 1 on hitachi systems])
AH_TEMPLATE([HPUX],         [Define to 1 on HP-UX systems])
AH_TEMPLATE([HPUX64],       [Define to 1 on 64bit HP-UX systems])
AH_TEMPLATE([IBM],          [Define to 1 on IBM SP systems])
AH_TEMPLATE([IBM64],        [Define to 1 on 64bit IBM SP systems])
AH_TEMPLATE([INTERIX],      [Define to 1 on ??? systems])
AH_TEMPLATE([LAPI],         [Define to 1 on IBM systems with LAPI])
AH_TEMPLATE([LAPI64],       [Define to 1 on 64bit IBM systems with LAPI])
AH_TEMPLATE([LINUX],        [Define to 1 on generic Linux systems])
AH_TEMPLATE([LINUX64],      [Define to 1 on generic 64bit Linux systems])
AH_TEMPLATE([MACX],         [Define to 1 on OSX systems])
AH_TEMPLATE([MACX64],       [Define to 1 on 64bit OSX systems])
AH_TEMPLATE([NEC],          [Define to 1 on NEC systems])
AH_TEMPLATE([NEC64],        [Define to 1 on 64bit NEC systems])
AH_TEMPLATE([SGI],          [Define to 1 on ??? systems])
AH_TEMPLATE([SGI_N32],      [Define to 1 on ??? systems])
AH_TEMPLATE([SGITFP],       [Define to 1 on ??? systems])
AH_TEMPLATE([SOLARIS],      [Define to 1 on Solaris systems])
AH_TEMPLATE([SOLARIS64],    [Define to 1 on 64bit Solaris systems])
AC_REQUIRE([AC_CANONICAL_BUILD])
AC_REQUIRE([AC_CANONICAL_HOST])
AC_CACHE_CHECK([for TARGET base (64bit-ness checked later)],
[ga_cv_target_base],
[ga_cv_target_base=UNKNOWN
AS_IF([test "x$ga_cv_target_base" = xUNKNOWN],
    [AS_IF([test -f /bgsys/drivers/ppcfloor/arch/include/common/bgp_personality.h],
        [ga_cv_target_base=BGP])])
AS_IF([test "x$ga_cv_target_base" = xUNKNOWN],
    [AS_IF([test -d /bgl/BlueLight/ppcfloor/bglsys/include],
        [ga_cv_target_base=BGL])])
AS_IF([test "x$ga_cv_target_base" = xUNKNOWN],
    [AS_CASE([$host],
        [*bgl*],            [ga_cv_target_base=BGL],
        [*bgp*],            [ga_cv_target_base=BGP],
        #[TODO],            [ga_cv_target_base=CATAMOUNT],
        #[TODO],            [ga_cv_target_base=CRAY_XT],
        [*cygwin*],         [ga_cv_target_base=CYGWIN],
        [*fujitsu*],        [ga_cv_target_base=FUJITSU_VPP],
        [*hpux*],           [ga_cv_target_base=HPUX],
        [*ibm*],            [ga_cv_target_base=IBM],
        #[TODO],            [ga_cv_target_base=LAPI],
        [*linux*],          [ga_cv_target_base=LINUX],
        [*darwin*],         [ga_cv_target_base=MACX],
        [*apple*],          [ga_cv_target_base=MACX],
        [*superux*],        [ga_cv_target_base=NEC],
        [*solaris*],        [ga_cv_target_base=SOLARIS])])
])dnl
AC_DEFINE_UNQUOTED([$ga_cv_target_base], [1],
    [define if this is the TARGET irregardless of whether it is 32/64 bits])
# A horrible hack that should go away somehow...
dnl # Only perform this hack for ARMCI build.
dnl AS_IF([test "x$ARMCI_TOP_BUILDDIR" != x], [
    AC_CACHE_CHECK([whether we think this system is what we call SYSV],
    [ga_cv_sysv],
    [AS_CASE([$ga_cv_target_base],
        [SUN|SOLARIS|SGI|SGI_N32|SGITFP|HPUX|IBM|DECOSF|LINUX|INTERIX|NEC|LAPI],
            [ga_cv_sysv=yes],
        [ga_cv_sysv=no])
    ])
    AS_IF([test x$ga_cv_sysv = xyes],
        [AC_DEFINE([SYSV], [1],
            [Define if we want this system to use SYSV shared memory])])
dnl ])
# Hopefully these will never be used and we can remove them soon.
AM_CONDITIONAL([BGL],          [test "$ga_cv_target_base" = BGL])
AM_CONDITIONAL([BGP],          [test "$ga_cv_target_base" = BGP])
AM_CONDITIONAL([CATAMOUNT],    [test "$ga_cv_target_base" = CATAMOUNT])
AM_CONDITIONAL([CRAY_SV1],     [test "$ga_cv_target_base" = CRAY_SV1])
AM_CONDITIONAL([CRAY_SV2],     [test "$ga_cv_target_base" = CRAY_SV2])
AM_CONDITIONAL([CRAY_T3E],     [test "$ga_cv_target_base" = CRAY_T3E])
AM_CONDITIONAL([CRAY_XT],      [test "$ga_cv_target_base" = CRAY_XT])
AM_CONDITIONAL([CRAY_YMP],     [test "$ga_cv_target_base" = CRAY_YMP])
AM_CONDITIONAL([CYGNUS],       [test "$ga_cv_target_base" = CYGNUS])
AM_CONDITIONAL([CYGWIN],       [test "$ga_cv_target_base" = CYGWIN])
AM_CONDITIONAL([DECOSF],       [test "$ga_cv_target_base" = DECOSF])
AM_CONDITIONAL([FUJITSU_VPP],  [test "$ga_cv_target_base" = FUJITSU_VPP])
AM_CONDITIONAL([HITACHI],      [test "$ga_cv_target_base" = HITACHI])
AM_CONDITIONAL([HPUX],         [test "$ga_cv_target_base" = HPUX])
AM_CONDITIONAL([IBM],          [test "$ga_cv_target_base" = IBM])
AM_CONDITIONAL([INTERIX],      [test "$ga_cv_target_base" = INTERIX])
AM_CONDITIONAL([LAPI],         [test "$ga_cv_target_base" = LAPI])
AM_CONDITIONAL([LINUX],        [test "$ga_cv_target_base" = LINUX])
AM_CONDITIONAL([MACX],         [test "$ga_cv_target_base" = MACX])
AM_CONDITIONAL([NEC],          [test "$ga_cv_target_base" = NEC])
AM_CONDITIONAL([SGI],          [test "$ga_cv_target_base" = SGI])
AM_CONDITIONAL([SGI_N32],      [test "$ga_cv_target_base" = SGI_N32])
AM_CONDITIONAL([SGITFP],       [test "$ga_cv_target_base" = SGITFP])
AM_CONDITIONAL([SOLARIS],      [test "$ga_cv_target_base" = SOLARIS])
])dnl


# GA_TARGET64()
# -------------
# Checking for 64bit platforms requires checking sizeof void*.
# That's easy, but doing it too soon causes AC_PROG_F77/C/CXX to get expanded
# too soon, and we want to expand those with a better list of compilers
# based on our current TARGET. Therefore, we must do this 64bit test later.
AC_DEFUN([GA_TARGET64],
[AC_REQUIRE([GA_TARGET])
AC_COMPUTE_INT([ga_target64_sizeof_voidp],
    [(long int) (sizeof (void*))])
AC_CACHE_CHECK([for TARGET 64bit-ness], [ga_cv_target],
[AS_IF([test x$ga_target64_sizeof_voidp = x8],
    [ga_cv_target=${ga_cv_target_base}64],
    [ga_cv_target=$ga_cv_target_base])])
AC_DEFINE_UNQUOTED([$ga_cv_target], [1],
    [define if this is the TARGET, 64bit-specific])
])dnl
