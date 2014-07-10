# COMEX_C_BLAS_TEST
# -----------------
# Generate C conftest for BLAS.
AC_DEFUN([COMEX_C_BLAS_TEST], [
AC_LANG_CONFTEST([AC_LANG_PROGRAM(
[#ifdef __cplusplus
extern "C" {
#endif
char $caxpy ();
char $daxpy ();
char $saxpy ();
char $zaxpy ();
char $ccopy ();
char $dcopy ();
char $scopy ();
char $zcopy ();
#ifdef __cplusplus
}
#endif
],
[[char caxpy_result = $caxpy ();
char daxpy_result = $daxpy ();
char saxpy_result = $saxpy ();
char zaxpy_result = $zaxpy ();
char ccopy_result = $ccopy ();
char dcopy_result = $dcopy ();
char scopy_result = $scopy ();
char zcopy_result = $zcopy ();
]])])
])

# COMEX_RUN_BLAS_TEST
# -------------------
# Test the linker.
# Clears BLAS_LIBS on failure.  Sets comex_blas_ok=yes on success.
AC_DEFUN([COMEX_RUN_BLAS_TEST], [
    AC_LANG_PUSH([C])
    comex_blas_ok=no
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=caxpy
         daxpy=daxpy
         saxpy=saxpy
         zaxpy=zaxpy
         ccopy=ccopy
         dcopy=dcopy
         scopy=scopy
         zcopy=zcopy
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=caxpy_
         daxpy=daxpy_
         saxpy=saxpy_
         zaxpy=zaxpy_
         ccopy=ccopy_
         dcopy=dcopy_
         scopy=scopy_
         zcopy=zcopy_
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=caxpy__
         daxpy=daxpy__
         saxpy=saxpy__
         zaxpy=zaxpy__
         ccopy=ccopy__
         dcopy=dcopy__
         scopy=scopy__
         zcopy=zcopy__
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=CGEMM
         daxpy=DGEMM
         saxpy=SGEMM
         zaxpy=ZGEMM
         ccopy=CCOPY
         dcopy=DCOPY
         scopy=SCOPY
         zcopy=ZCOPY
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=CGEMM_
         daxpy=DGEMM_
         saxpy=SGEMM_
         zaxpy=ZGEMM_
         ccopy=CCOPY_
         dcopy=DCOPY_
         scopy=SCOPY_
         zcopy=ZCOPY_
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=CGEMM__
         daxpy=DGEMM__
         saxpy=SGEMM__
         zaxpy=ZGEMM__
         ccopy=CCOPY__
         dcopy=DCOPY__
         scopy=SCOPY__
         zcopy=ZCOPY__
         COMEX_C_BLAS_TEST()
         AC_LINK_IFELSE([], [comex_blas_ok=yes], [BLAS_LIBS=])])
    AS_IF([test "x$comex_blas_ok" = xno],
        [caxpy=NOTFOUND
         daxpy=NOTFOUND
         saxpy=NOTFOUND
         zaxpy=NOTFOUND
         ccopy=NOTFOUND
         dcopy=NOTFOUND
         scopy=NOTFOUND
         zcopy=NOTFOUND])
    AC_LANG_POP([C])
])dnl

# COMEX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# -------------------------------------------------
# Parse the --with-blas argument and test for all *axpy routines. We use
# different tests depending on whether Fortran sources are enabled. There are
# many flavors of BLAS that we test for explicitly, although the list could
# probably be reduced based on currently available systems.
AC_DEFUN([COMEX_BLAS], [
blas_size=4
blas_size_hack=no
AC_ARG_WITH([blas],
    [AS_HELP_STRING([--with-blas[[=ARG]]],
        [use external BLAS library; attempt to detect sizeof(INTEGER)])],
    [blas_size_hack=yes])
AC_ARG_WITH([blas4],
    [AS_HELP_STRING([--with-blas4[[=ARG]]],
        [use external BLAS library compiled with sizeof(INTEGER)==4])],
    [blas_size=4; with_blas="$with_blas4"])
AC_ARG_WITH([blas8],
    [AS_HELP_STRING([--with-blas8[[=ARG]]],
        [use external BLAS library compiled with sizeof(INTEGER)==8])],
    [blas_size=8; with_blas="$with_blas8"])

comex_blas_ok=no
AS_IF([test "x$with_blas" = xno], [comex_blas_ok=skip])

# Parse --with-blas argument. Clear previous values first.
BLAS_LIBS=
BLAS_LDFLAGS=
BLAS_CPPFLAGS=
COMEX_ARG_PARSE([with_blas], [BLAS_LIBS], [BLAS_LDFLAGS], [BLAS_CPPFLAGS])

comex_save_LIBS="$LIBS"
comex_save_LDFLAGS="$LDFLAGS";     LDFLAGS="$BLAS_LDFLAGS $LDFLAGS"
comex_save_CPPFLAGS="$CPPFLAGS";   CPPFLAGS="$BLAS_CPPFLAGS $CPPFLAGS"

AC_MSG_NOTICE([Attempting to locate BLAS library])

# First, check environment/command-line variables.
# If failed, erase BLAS_LIBS but maintain BLAS_LDFLAGS and BLAS_CPPFLAGS.
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS with user-supplied flags])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AS_IF([test $comex_blas_ok = yes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*ilp64*],  [blas_size=8],     # Intel MKL
                [*_int64*], [blas_size=8])])]) # AMD ACML
     AC_MSG_RESULT([$comex_blas_ok])])

# AMD Core Math Library (ACML)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in AMD Core Math Library])
     # add -lacml to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*acml*], [], [BLAS_LIBS="-lacml"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AS_IF([test "x$comex_blas_ok" = xyes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*_int64*], [blas_size=8])])])
     AC_MSG_RESULT([$comex_blas_ok])])

# Intel MKL library
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Intel Math Kernel Library])
     # add -lmkl to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*mkl*], [], [BLAS_LIBS="-lmkl"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AS_IF([test "x$comex_blas_ok" = xyes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*ilp64*], [blas_size=8])])])
     AC_MSG_RESULT([$comex_blas_ok])])

# ATLAS library (http://math-atlas.sourceforge.net/)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in ATLAS])
     AS_IF([test "x$enable_f77" = xno],
        [# add -lcblas if needed but missing from LIBS
         AS_CASE([$LIBS], [*cblas*], [], [BLAS_LIBS="-lcblas"])],
        [# add -lf77blas if needed but missing from LIBS
         AS_CASE([$LIBS], [*f77blas*], [], [BLAS_LIBS="-lf77blas"])])
     # add -latlas if needed but missing from LIBS
     AS_CASE([$LIBS], [*atlas*], [], [BLAS_LIBS="$BLAS_LIBS -latlas"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# PhiPACK libraries (requires generic BLAS lib, too)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in PhiPACK libraries])
     # add -lblas to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*blas*], [], [BLAS_LIBS="-lblas"])
     # add -ldgemm to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*dgemm*], [], [BLAS_LIBS="-ldgemm $BLAS_LIBS"])
     # add -lsgemm to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*sgemm*], [], [BLAS_LIBS="-lsgemm $BLAS_LIBS"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# Apple Accelerate.framework
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Apple Accelerate.framework])
     # add -framework Accelerate to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*Accelerate*], [], [BLAS_LIBS="-framework Accelerate"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# Apple vecLib.framework
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Apple vecLib.framework])
     # add -framework vecLib to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*vecLib*], [], [BLAS_LIBS="-framework vecLib"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# Alpha CXML library (CXML stands for Compaq Extended Math Library)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Alpha CXML library])
     # add -lcxml to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*cxml*], [], [BLAS_LIBS="-lcxml"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AS_IF([test $comex_blas_ok = no],
        [# add -lcxml to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*cxml*], [], [BLAS_LIBS="-lcxml"])
         # add -lcpml to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*cpml*], [], [BLAS_LIBS="$BLAS_LIBS -lcpml"])
         LIBS="$BLAS_LIBS $LIBS"
         COMEX_RUN_BLAS_TEST()
         LIBS="$comex_save_LIBS"])
     AC_MSG_RESULT([$comex_blas_ok])])

# Alpha DXML library (now called CXML, see above)

# Sun Performance library
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Sun Performance Library])
     # add -xlic_lib=sunperf to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*sunperf*], [], [BLAS_LIBS="-xlic_lib=sunperf"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AS_IF([test $comex_blas_ok = no],
        [# add -xlic_lib=sunperf to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*sunperf*], [], [BLAS_LIBS="-xlic_lib=sunperf"])
         # add -lsunmath to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*sunmath*], [], [BLAS_LIBS="$BLAS_LIBS -lsunmath"])
         LIBS="$BLAS_LIBS $LIBS"
         COMEX_RUN_BLAS_TEST()
         LIBS="$comex_save_LIBS"])
     AC_MSG_RESULT([$comex_blas_ok])])

# SCSL library (SCSL stands for SGI/Cray Scientific Library)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in SGI/Cray Scientific Library])
     # add -lscs to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*scs*], [], [BLAS_LIBS="-lscs"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# SGIMATH library
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in SGIMATH library])
     # add -lcomplib.sgimath to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*complib.sgimath*], [], [BLAS_LIBS="-lcomplib.sgimath"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# IBM ESSL library (requires generic BLAS lib, too)
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in IBM ESSL library])
     # add -lessl to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*essl*], [], [BLAS_LIBS="-lessl"])
     # add -lblas to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*blas*], [], [BLAS_LIBS="$BLAS_LIBS -lblas"])
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

# Generic BLAS library
AS_IF([test $comex_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in generic library])
     BLAS_LIBS="-lblas"
     LIBS="$BLAS_LIBS $LIBS"
     COMEX_RUN_BLAS_TEST()
     LIBS="$comex_save_LIBS"
     AC_MSG_RESULT([$comex_blas_ok])])

CPPFLAGS="$comex_save_CPPFLAGS"
LDFLAGS="$comex_save_LDFLAGS"

AC_SUBST([BLAS_LIBS])
AC_SUBST([BLAS_LDFLAGS])
AC_SUBST([BLAS_CPPFLAGS])

# Tests are complete. Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $comex_blas_ok = yes],
    [have_blas=1
     $1],
    [AC_MSG_WARN([BLAS library not found, using internal BLAS])
     blas_size=4
     have_blas=0
     $2])
AC_DEFINE_UNQUOTED([HAVE_BLAS], [$have_blas],
    [Define to 1 if using external BLAS library])
AC_DEFINE_UNQUOTED([BLAS_SIZE], [$blas_size],
    [Define to sizeof(INTEGER) used to compile BLAS])
AC_DEFINE_UNQUOTED([BLAS_CAXPY], [$caxpy],
    [Define to name of caxpy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_DAXPY], [$daxpy],
    [Define to name of daxpy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_SAXPY], [$saxpy],
    [Define to name of saxpy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_ZAXPY], [$zaxpy],
    [Define to name of zaxpy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_CCOPY], [$ccopy],
    [Define to name of ccopy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_DCOPY], [$dcopy],
    [Define to name of dcopy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_SCOPY], [$scopy],
    [Define to name of scopy routine to call from C])
AC_DEFINE_UNQUOTED([BLAS_ZCOPY], [$zcopy],
    [Define to name of zcopy routine to call from C])
AM_CONDITIONAL([HAVE_BLAS], [test $comex_blas_ok = yes])
])dnl COMEX_BLAS
