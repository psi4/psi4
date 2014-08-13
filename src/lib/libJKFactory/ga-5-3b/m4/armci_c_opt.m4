# ARMCI_C_OPT()
# -------------
# Determine TARGET-/compiler-specific CFLAGS for optimization.
AC_DEFUN([ARMCI_C_OPT], [
AC_REQUIRE([GA_TARGET64])
AC_REQUIRE([GA_ENABLE_OPT])
AC_REQUIRE([GA_ARMCI_NETWORK])
AC_ARG_VAR([ARMCI_COPT], [ARMCI C optimization flags])
AC_CACHE_CHECK([for specific C optimizations], [armci_cv_c_opt], [
AS_IF([test "x$ARMCI_COPT" != x], [armci_cv_c_opt="$ARMCI_COPT"], [armci_cv_c_opt=])
AS_IF([test "x$armci_cv_c_opt" = x && test "x$enable_opt" = xyes], [
AS_CASE([$ga_cv_target:$ga_cv_c_compiler_vendor:$host_cpu:$ga_armci_network],
[BGL:*:*:*],                [armci_cv_c_opt="-O0"],
[BGP:ibm:*:*],              [armci_cv_c_opt="-O3 -qstrict -qarch=450 -qtune=450"],
[BGP:gnu:*:*],              [armci_cv_c_opt="-O2"],
[CATAMOUNT:*:*:*],          [armci_cv_c_opt=],
[CRAY_XT:*:*:*],            [armci_cv_c_opt=],
[CYGWIN:*:*:*],             [armci_cv_c_opt="-malign-double"],
[FUJITSU_VPP64:*:*:*],      [armci_cv_c_opt="-x100"],
[FUJITSU_VPP:*:*:*],        [armci_cv_c_opt="-x100 -KA32"],
[HPUX64:*:*:*],             [armci_cv_c_opt="-Ae"],
[HPUX64:*:ia64:*],          [armci_cv_c_opt="-Ae"],
[HPUX:*:*:*],               [armci_cv_c_opt="-Ae"],
[IBM64:*:*:*],              [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[IBM:*:*:*],                [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LAPI64:*:*:*],             [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LAPI:*:*:*],               [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LINUX64:fujitsu:ia64:*],   [armci_cv_c_opt="-Kfast"],
[LINUX64:fujitsu:x86_64:*], [armci_cv_c_opt="-Kfast"],
[LINUX64:gnu:ia64:*],       [armci_cv_c_opt="-O0 -g"],
[LINUX64:gnu:x86_64:*],     [armci_cv_c_opt="-O3 -funroll-loops"],
[LINUX64:ibm:powerpc64:*],  [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LINUX64:ibm:ppc64:*],      [armci_cv_c_opt="-O3 -qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LINUX64:ibm:x86_64:*],     [armci_cv_c_opt=""],
[LINUX64:intel:ia64:*],     [armci_cv_c_opt="-w1"],
[LINUX64:unknown:alpha:*],  [armci_cv_c_opt="-assume no2underscore -fpe3 -check nooverflow -assume accuracy_sensitive -check nopower -check nounderflow"],
[LINUX:fujitsu:*:*],        [armci_cv_c_opt="-Kfast"],
[LINUX:gnu:686:*],          [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro -malign-double"],
[LINUX:gnu:686:MELLANOX],   [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro"],
[LINUX:gnu:686:OPENIB],     [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro"],
[LINUX:gnu:786:*],          [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro -malign-double"],
[LINUX:gnu:786:MELLANOX],   [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro"],
[LINUX:gnu:786:OPENIB],     [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -march=pentiumpro"],
[LINUX:gnu:x86:*],          [armci_cv_c_opt="-O2 -finline-functions -funroll-loops -malign-double"],
[LINUX:gnu:x86:MELLANOX],   [armci_cv_c_opt="-O2 -finline-functions -funroll-loops "],
[LINUX:gnu:x86:OPENIB],     [armci_cv_c_opt="-O2 -finline-functions -funroll-loops "],
[LINUX:ibm:*:*],            [armci_cv_c_opt="-q32"],
[LINUX:intel:*:*],          [armci_cv_c_opt="-O3 -prefetch"],
[MACX64:*:*:*],             [armci_cv_c_opt=],
[MACX:*:*:*],               [armci_cv_c_opt=],
[NEC64:*:*:*],              [armci_cv_c_opt="-Cvsafe -size_t64"],
[NEC:*:*:*],                [armci_cv_c_opt="-Cvsafe"],
[SOLARIS64:fujitsu:*:*],    [armci_cv_c_opt="-Kfast -KV9FMADD -x0"],
[SOLARIS64:gnu:*:*],        [armci_cv_c_opt="-dalign"],
[SOLARIS64:gnu:i386:*],     [armci_cv_c_opt="-dalign -xarch=amd64"],
[SOLARIS:fujitsu:*:*],      [armci_cv_c_opt="-Kfast -KV8PFMADD -x0"],
[SOLARIS:gnu:*:*],          [armci_cv_c_opt="-dalign"],
[SOLARIS:gnu:i386:*],       [armci_cv_c_opt="-dalign -xarch=sse2"],
                            [armci_cv_c_opt=])
])])
AC_SUBST([ARMCI_COPT],  [$armci_cv_c_opt])
])dnl
