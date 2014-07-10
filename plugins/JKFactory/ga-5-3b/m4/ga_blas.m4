# GA_F77_BLAS_TEST
# ----------------
# Generate Fortran 77 conftest for BLAS.
AC_DEFUN([GA_F77_BLAS_TEST], [AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
[[      implicit none
      INTEGER M, N, K, LDA, LDB, LDC
      COMPLEX CA(20,40), CB(20,30), CC(40,30), Calpha, Cbeta
      DOUBLE COMPLEX ZA(20,40), ZB(20,30), ZC(40,30), Zalpha, Zbeta
      REAL SA(20,40), SB(20,30), SC(40,30), Salpha, Sbeta
      DOUBLE PRECISION DA(20,40), DB(20,30), DC(40,30), Dalpha, Dbeta
      external CGEMM
      external ZGEMM
      external SGEMM
      external DGEMM
      M = 10
      N = 20
      K = 15
      LDA = 20
      LDB = 20
      LDC = 40
      Calpha = 2.0
      Cbeta = 2.0
      Zalpha = 2.0
      Zbeta = 2.0
      Salpha = 2.0
      Sbeta = 2.0
      Dalpha = 2.0
      Dbeta = 2.0
      CALL CGEMM ('T','N',M,N,K,Calpha,CA,LDA,CB,LDB,Cbeta,CC,LDC)
      CALL ZGEMM ('T','N',M,N,K,Zalpha,ZA,LDA,ZB,LDB,Zbeta,ZC,LDC)
      CALL SGEMM ('T','N',M,N,K,Salpha,SA,LDA,SB,LDB,Sbeta,SC,LDC)
      CALL DGEMM ('T','N',M,N,K,Dalpha,DA,LDA,DB,LDB,Dbeta,DC,LDC)]])])
])

# GA_C_BLAS_TEST
# --------------
# Generate C conftest for BLAS.
AC_DEFUN([GA_C_BLAS_TEST], [AC_LANG_CONFTEST([AC_LANG_PROGRAM(
[#ifdef __cplusplus
extern "C" {
#endif
char cgemm ();
char dgemm ();
char sgemm ();
char zgemm ();
#ifdef __cplusplus
}
#endif
],
[[char cresult =  cgemm ();
char dresult =  dgemm ();
char sresult =  sgemm ();
char zresult =  zgemm ();
]])])
])

# GA_C_RUN_BLAS_TEST
# ------------------
# Test the C linker.
# Clears BLAS_LIBS on failure.  Sets ga_blas_ok=yes on success.
AC_DEFUN([GA_C_RUN_BLAS_TEST], [
   AC_LANG_PUSH([C])
   GA_C_BLAS_TEST()
   AC_LINK_IFELSE([], [ga_blas_ok=yes], [BLAS_LIBS=])
   AC_LANG_POP([C])
])dnl

# GA_F77_RUN_BLAS_TEST
# --------------------
# Test the Fortran 77 linker.
# Clears BLAS_LIBS on failure.  Sets ga_blas_ok=yes on success.
AC_DEFUN([GA_F77_RUN_BLAS_TEST], [
   AC_LANG_PUSH([Fortran 77])
   GA_F77_BLAS_TEST()
   AC_LINK_IFELSE([], [ga_blas_ok=yes], [BLAS_LIBS=])
   AC_LANG_POP([Fortran 77])
])dnl

# GA_RUN_BLAS_TEST
# ----------------
# Test the linker.
# Clears BLAS_LIBS on failure.  Sets ga_blas_ok=yes on success.
AC_DEFUN([GA_RUN_BLAS_TEST], [
AS_IF([test "x$enable_f77" = xno],
   [AC_LANG_PUSH([C])
    GA_C_BLAS_TEST()
    AC_LINK_IFELSE([], [ga_blas_ok=yes], [BLAS_LIBS=])
    AC_LANG_POP([C])],
   [AC_LANG_PUSH([Fortran 77])
    GA_F77_BLAS_TEST()
    AC_LINK_IFELSE([], [ga_blas_ok=yes], [BLAS_LIBS=])
    AC_LANG_POP([Fortran 77])])
])dnl

# GA_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# -------------------------------------------------
# Parse the --with-blas argument and test for all *gemm routines. We use
# different tests depending on whether Fortran sources are enabled. There are
# many flavors of BLAS that we test for explicitly, although the list could
# probably be reduced based on currently available systems.
#
# Apparently certain compilers on BGP define sgemm and dgemm, so we must
# test for a different BLAS routine. cgemm seems okay.
AC_DEFUN([GA_BLAS],
[AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
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

ga_blas_ok=no
AS_IF([test "x$with_blas" = xno], [ga_blas_ok=skip])

# Parse --with-blas argument. Clear previous values first.
BLAS_LIBS=
BLAS_LDFLAGS=
BLAS_CPPFLAGS=
GA_ARG_PARSE([with_blas], [BLAS_LIBS], [BLAS_LDFLAGS], [BLAS_CPPFLAGS])

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(cgemm)
AC_F77_FUNC(dgemm)
AC_F77_FUNC(sgemm)
AC_F77_FUNC(zgemm)

ga_save_LIBS="$LIBS"
ga_save_LDFLAGS="$LDFLAGS";     LDFLAGS="$BLAS_LDFLAGS $LDFLAGS"
ga_save_CPPFLAGS="$CPPFLAGS";   CPPFLAGS="$BLAS_CPPFLAGS $CPPFLAGS"

AC_MSG_NOTICE([Attempting to locate BLAS library])

# First, check environment/command-line variables.
# If failed, erase BLAS_LIBS but maintain BLAS_LDFLAGS and BLAS_CPPFLAGS.
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS with user-supplied flags])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AS_IF([test $ga_blas_ok = yes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*ilp64*],  [blas_size=8],     # Intel MKL
                [*_int64*], [blas_size=8])])]) # AMD ACML
     AC_MSG_RESULT([$ga_blas_ok])])

# AMD Core Math Library (ACML)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in AMD Core Math Library])
     # add -lacml to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*acml*], [], [BLAS_LIBS="-lacml"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AS_IF([test "x$ga_blas_ok" = xyes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*_int64*], [blas_size=8])])])
     AC_MSG_RESULT([$ga_blas_ok])])

# Intel MKL library
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Intel Math Kernel Library])
     # add -lmkl to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*mkl*], [], [BLAS_LIBS="-lmkl"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AS_IF([test "x$ga_blas_ok" = xyes],
        [AS_IF([test $blas_size_hack = yes],
            [AS_CASE(["$BLAS_LIBS $LIBS $LDFLAGS $CPPFLAGS"],
                [*ilp64*], [blas_size=8])])])
     AC_MSG_RESULT([$ga_blas_ok])])

# ATLAS library (http://math-atlas.sourceforge.net/)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in ATLAS])
     AS_IF([test "x$enable_f77" = xno],
        [# add -lcblas if needed but missing from LIBS
         AS_CASE([$LIBS], [*cblas*], [], [BLAS_LIBS="-lcblas"])],
        [# add -lf77blas if needed but missing from LIBS
         AS_CASE([$LIBS], [*f77blas*], [], [BLAS_LIBS="-lf77blas"])])
     # add -latlas if needed but missing from LIBS
     AS_CASE([$LIBS], [*atlas*], [], [BLAS_LIBS="$BLAS_LIBS -latlas"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# PhiPACK libraries (requires generic BLAS lib, too)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in PhiPACK libraries])
     # add -lblas to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*blas*], [], [BLAS_LIBS="-lblas"])
     # add -ldgemm to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*dgemm*], [], [BLAS_LIBS="-ldgemm $BLAS_LIBS"])
     # add -lsgemm to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*sgemm*], [], [BLAS_LIBS="-lsgemm $BLAS_LIBS"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# Apple Accelerate.framework
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Apple Accelerate.framework])
     # add -framework Accelerate to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*Accelerate*], [], [BLAS_LIBS="-framework Accelerate"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# Apple vecLib.framework
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Apple vecLib.framework])
     # add -framework vecLib to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*vecLib*], [], [BLAS_LIBS="-framework vecLib"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# Alpha CXML library (CXML stands for Compaq Extended Math Library)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Alpha CXML library])
     # add -lcxml to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*cxml*], [], [BLAS_LIBS="-lcxml"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AS_IF([test $ga_blas_ok = no],
        [# add -lcxml to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*cxml*], [], [BLAS_LIBS="-lcxml"])
         # add -lcpml to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*cpml*], [], [BLAS_LIBS="$BLAS_LIBS -lcpml"])
         LIBS="$BLAS_LIBS $LIBS"
         GA_RUN_BLAS_TEST()
         LIBS="$ga_save_LIBS"])
     AC_MSG_RESULT([$ga_blas_ok])])

# Alpha DXML library (now called CXML, see above)

# Sun Performance library
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in Sun Performance Library])
     # add -xlic_lib=sunperf to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*sunperf*], [], [BLAS_LIBS="-xlic_lib=sunperf"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AS_IF([test $ga_blas_ok = no],
        [# add -xlic_lib=sunperf to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*sunperf*], [], [BLAS_LIBS="-xlic_lib=sunperf"])
         # add -lsunmath to BLAS_LIBS if missing from LIBS
         AS_CASE([$LIBS], [*sunmath*], [], [BLAS_LIBS="$BLAS_LIBS -lsunmath"])
         LIBS="$BLAS_LIBS $LIBS"
         GA_RUN_BLAS_TEST()
         LIBS="$ga_save_LIBS"])
     AC_MSG_RESULT([$ga_blas_ok])])

# SCSL library (SCSL stands for SGI/Cray Scientific Library)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in SGI/Cray Scientific Library])
     # add -lscs to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*scs*], [], [BLAS_LIBS="-lscs"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# SGIMATH library
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in SGIMATH library])
     # add -lcomplib.sgimath to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*complib.sgimath*], [], [BLAS_LIBS="-lcomplib.sgimath"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# IBM ESSL library (requires generic BLAS lib, too)
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in IBM ESSL library])
     # add -lessl to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*essl*], [], [BLAS_LIBS="-lessl"])
     # add -lblas to BLAS_LIBS if missing from LIBS
     AS_CASE([$LIBS], [*blas*], [], [BLAS_LIBS="$BLAS_LIBS -lblas"])
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

# Generic BLAS library
AS_IF([test $ga_blas_ok = no],
    [AC_MSG_CHECKING([for BLAS in generic library])
     BLAS_LIBS="-lblas"
     LIBS="$BLAS_LIBS $LIBS"
     GA_RUN_BLAS_TEST()
     LIBS="$ga_save_LIBS"
     AC_MSG_RESULT([$ga_blas_ok])])

CPPFLAGS="$ga_save_CPPFLAGS"
LDFLAGS="$ga_save_LDFLAGS"

AC_SUBST([BLAS_LIBS])
AC_SUBST([BLAS_LDFLAGS])
AC_SUBST([BLAS_CPPFLAGS])

# Tests are complete. Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $ga_blas_ok = yes],
    [have_blas=1
     $1],
    [AC_MSG_WARN([BLAS library not found, using internal BLAS])
     blas_size=$ga_cv_f77_integer_size # reset blas integer size to desired
     have_blas=0
     $2])
AC_DEFINE_UNQUOTED([HAVE_BLAS], [$have_blas],
    [Define to 1 if using external BLAS library])
AC_DEFINE_UNQUOTED([BLAS_SIZE], [$blas_size],
    [Define to sizeof(INTEGER) used to compile BLAS])
AM_CONDITIONAL([HAVE_BLAS], [test $ga_blas_ok = yes])
])dnl GA_BLAS
