# GA_C_OPT()
# ----------
# Determine TARGET-/compiler-specific CFLAGS for optimization.
AC_DEFUN([GA_C_OPT], [
AC_REQUIRE([GA_TARGET64])
AC_REQUIRE([GA_ENABLE_OPT])
AC_ARG_VAR([GA_COPT], [GA C optimization flags])
AC_CACHE_CHECK([for specific C optimizations], [ga_cv_c_opt], [
AS_IF([test "x$GA_COPT" != x], [ga_cv_c_opt="$GA_COPT"], [ga_cv_c_opt=])
AS_IF([test "x$ga_cv_c_opt" = x && test "x$enable_opt" = xyes], [
AS_CASE([$ga_cv_target:$ga_cv_c_compiler_vendor:$host_cpu],
[BGL:*:*],                  [ga_cv_c_opt="-O0"],
[BGP:ibm:*],                [ga_cv_c_opt="-O3 -qstrict -qarch=450 -qtune=450"],
[BGP:gnu:*],                [ga_cv_c_opt="-O2"],
[CATAMOUNT:*:*],            [ga_cv_c_opt=],
[CRAY_XT:*:*],              [ga_cv_c_opt=],
[CYGWIN:*:*],               [ga_cv_c_opt=],
[FUJITSU_VPP64:*:*],        [ga_cv_c_opt=],
[FUJITSU_VPP:*:*],          [ga_cv_c_opt="-KA32"],
[HPUX64:*:*],               [ga_cv_c_opt="-Ae"],
[HPUX64:*:ia64],            [ga_cv_c_opt="-Ae"],
[HPUX:*:*],                 [ga_cv_c_opt="-Ae"],
[IBM64:*:*],                [ga_cv_c_opt=],
[IBM:*:*],                  [ga_cv_c_opt=],
[LAPI64:*:*],               [ga_cv_c_opt=],
[LAPI:*:*],                 [ga_cv_c_opt=],
[LINUX64:fujitsu:ia64],     [ga_cv_c_opt="-Kfast"],
[LINUX64:fujitsu:x86_64],   [ga_cv_c_opt="-Kfast"],
[LINUX64:gnu:ia64],         [ga_cv_c_opt="-O3 -funroll-loops"],
[LINUX64:gnu:powerpc64],    [ga_cv_c_opt="-funroll-loops"],
[LINUX64:gnu:ppc64],        [ga_cv_c_opt="-funroll-loops"],
[LINUX64:gnu:x86_64],       [ga_cv_c_opt="-O2 -funroll-loops"],
[LINUX64:ibm:powerpc64],    [ga_cv_c_opt="-qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LINUX64:ibm:ppc64],        [ga_cv_c_opt="-qinline=100 -qstrict -qarch=auto -qtune=auto"],
[LINUX64:ibm:x86_64],       [ga_cv_c_opt=],
[LINUX64:intel:ia64],       [ga_cv_c_opt="-fno-alias -ftz"],
[LINUX:fujitsu:*],          [ga_cv_c_opt="-Kfast"],
[LINUX:gnu:786],            [ga_cv_c_opt="-O2 -funroll-loops -malign-double"],
[LINUX:gnu:*],              [ga_cv_c_opt="-O2 -funroll-loops"],
[LINUX:gnu:x86],            [ga_cv_c_opt="-O2 -funroll-loops -malign-double"],
[LINUX:ibm:*],              [ga_cv_c_opt="-q32"],
[LINUX:intel:*],            [ga_cv_c_opt="-O3 -prefetch"],
[MACX64:*:*],               [ga_cv_c_opt=],
[MACX:*:*],                 [ga_cv_c_opt=],
[NEC64:*:*],                [ga_cv_c_opt="-Cvsafe -size_t64"],
[NEC:*:*],                  [ga_cv_c_opt="-Cvsafe"],
[SOLARIS64:fujitsu:*],      [ga_cv_c_opt="-Kfast -KV9FMADD"],
[SOLARIS64:gnu:*],          [ga_cv_c_opt="-dalign -xarch=v9"],
[SOLARIS64:gnu:i386],       [ga_cv_c_opt="-dalign -xarch=amd64"],
[SOLARIS:fujitsu:*],        [ga_cv_c_opt="-Kfast -KV8PFMADD"],
[SOLARIS:gnu:*],            [ga_cv_c_opt="-dalign"],
[SOLARIS:gnu:i386],         [ga_cv_c_opt="-dalign -xarch=sse2"],
                            [ga_cv_c_opt=])
])])
AC_SUBST([GA_COPT], [$ga_cv_c_opt])
])dnl
