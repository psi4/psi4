# COMEX_PROG_MPICC
# ----------------
# If desired, replace CC with MPICC while searching for a C compiler.
#
# Known C compilers
#  cc       generic compiler name
#  ccc      Fujitsu ?? old Cray ??
#  cl
#  ecc      Intel on IA64 ??
#  gcc      GNU
#  icc      Intel
#  bgxlc    Intel on BG/P
#  bgxlc_r  Intel on BG/P, thread safe
#  xlc      Intel
#  xlc_r    Intel, thread safe
#  pgcc     Portland Group
#  pathcc   PathScale
#  sxcc     NEC SX
#  fcc      Fujitsu
#  opencc   AMD's x86 open64
#  suncc    Sun's Studio
#  craycc   Cray
#
# Known MPI C compilers:
#  mpicc
#  mpixlc_r
#  mpixlc
#  hcc
#  mpxlc_r
#  mpxlc
#  sxmpicc  NEC SX
#  mpifcc   Fujitsu
#  mpgcc
#  mpcc
#  cmpicc
#  cc
#
AC_DEFUN([COMEX_PROG_MPICC],
[AC_ARG_VAR([MPICC], [MPI C compiler])
AS_CASE([$comex_cv_target_base],
[BGP],  [comex_mpicc_pref=mpixlc_r; comex_cc_pref=bgxlc_r],
[])
# In the case of using MPI wrappers, set CC=MPICC since CC will override
# absolutely everything in our list of compilers.
# Save CC, just in case.
AS_IF([test x$with_mpi_wrappers = xyes],
    [AS_IF([test "x$CC" != "x$MPICC"], [comex_orig_CC="$CC"])
     AS_CASE([x$CC:x$MPICC],
        [x:x],  [],
        [x:x*], [CC="$MPICC"],
        [x*:x],
[AC_MSG_WARN([MPI compilers desired but CC is set while MPICC is unset.])
 AC_MSG_WARN([CC will be ignored during compiler selection, but will be])
 AC_MSG_WARN([tested first during MPI compiler unwrapping. Perhaps you])
 AC_MSG_WARN([meant to set MPICC instead of or in addition to CC?])
 CC=],
        [x*:x*], 
[AS_IF([test "x$CC" != "x$MPICC"],
[AC_MSG_WARN([MPI compilers desired, MPICC and CC are set, and MPICC!=CC.])
 AC_MSG_WARN([Choosing MPICC over CC.])
 AC_MSG_WARN([CC will be tested first during MPI compiler unwrapping.])])
 CC="$MPICC"],
[AC_MSG_ERROR([CC/MPICC case failure])])])
comex_cc="bgxlc_r bgxlc xlc_r xlc pgcc pathcc icc sxcc fcc opencc suncc craycc gcc cc ecc cl ccc"
comex_mpicc="mpicc mpixlc_r mpixlc hcc mpxlc_r mpxlc sxmpicc mpifcc mpgcc mpcc cmpicc cc"
AS_IF([test x$with_mpi_wrappers = xyes],
    [CC_TO_TEST="$comex_mpicc_pref $comex_mpicc"],
    [CC_TO_TEST="$comex_cc_pref $comex_cc"])
AC_PROG_CC([$CC_TO_TEST])
])dnl
