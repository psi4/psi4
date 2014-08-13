# GA_F77_OPT()
# ------------
# Determine TARGET-/compiler-specific FFLAGS for optimization.
AC_DEFUN([GA_F77_OPT], [
AC_REQUIRE([GA_TARGET64])
AC_REQUIRE([GA_ENABLE_OPT])
AC_ARG_VAR([GA_FOPT], [GA Fortran 77 optimization flags])
AC_CACHE_CHECK([for specific Fortran optimizations], [ga_cv_f77_opt], [
AS_IF([test "x$GA_FOPT" != x], [ga_cv_f77_opt="$GA_FOPT"], [ga_cv_f77_opt=])
AS_IF([test "x$ga_cv_f77_opt" = x && test "x$enable_opt" = xyes], [
AS_CASE([$ga_cv_target:$ga_cv_f77_compiler_vendor:$host_cpu],
[BGL:*:*],                  [ga_cv_f77_opt="-O0"],
[BGP:ibm:*],                [ga_cv_f77_opt="-O3 -qstrict -qarch=450 -qtune=450"],
[BGP:gnu:*],                [ga_cv_f77_opt="-O2"],
[CATAMOUNT:*:*],            [ga_cv_f77_opt="-O3"],
[CRAY_XT:*:*],              [ga_cv_f77_opt=],
[CYGWIN:*:*],               [ga_cv_f77_opt=],
[FUJITSU_VPP64:*:*],        [ga_cv_f77_opt="-Sw"],
[FUJITSU_VPP:*:*],          [ga_cv_f77_opt="-Sw -KA32"],
[HPUX64:*:*],               [ga_cv_f77_opt="-O1"],    
[HPUX64:*:ia64],            [ga_cv_f77_opt="-O1"],    
[HPUX:*:*],                 [ga_cv_f77_opt="-O1"],
[IBM64:*:*],                [ga_cv_f77_opt="-qarch=auto"],
[IBM:*:*],                  [ga_cv_f77_opt="-qarch=auto"],
[LAPI64:*:*],               [ga_cv_f77_opt="-qarch=auto"],
[LAPI:*:*],                 [ga_cv_f77_opt="-qarch=auto"],
[LINUX64:*:alpha],          [ga_cv_f77_opt="-align_dcommons -fpe3 -check nooverflow -assume accuracy_sensitive -check nopower -check nounderflow"],
[LINUX64:fujitsu:ia64],     [ga_cv_f77_opt="-Kfast -X9 -Am -fw"],
[LINUX64:fujitsu:x86_64],   [ga_cv_f77_opt="-Kfast -X9 -Am -fw"],
[LINUX64:gnu:x86_64],       [ga_cv_f77_opt="-O"],
[LINUX64:ibm:x86_64],       [ga_cv_f77_opt=],
[LINUX64:intel:ia64],       [ga_cv_f77_opt="-cm -w90 -w95 -align"],
[LINUX64:intel:powerpc64],  [ga_cv_f77_opt=],
[LINUX64:intel:ppc64],      [ga_cv_f77_opt=],
[LINUX64:intel:x86_64],     [ga_cv_f77_opt="-O3 -w -cm -xW -tpp7"],
[LINUX64:portland:x86_64],  [ga_cv_f77_opt="-Mdalign"],
[LINUX:fujitsu:*],          [ga_cv_f77_opt="-Kfast -X9 -Am -fw"],
[LINUX:gnu:786],            [ga_cv_f77_opt="-O2 -funroll-loops -malign-double"],
[LINUX:gnu:*],              [ga_cv_f77_opt="-O2 -funroll-loops"],
[LINUX:gnu:x86],            [ga_cv_f77_opt="-O2 -funroll-loops -malign-double"],
[LINUX:ibm:*],              [ga_cv_f77_opt="-q32"],
[LINUX:intel:*],            [ga_cv_f77_opt="-O3 -prefetch -w -cm"],
[LINUX:portland:*],         [ga_cv_f77_opt="-Mdalign -Minform,warn -Mnolist -Minfo=loop -Munixlogical"],
[MACX64:ibm:*],             [ga_cv_f77_opt=],
[MACX64:intel:*],           [ga_cv_f77_opt="-O3 -prefetch -w -cm"],
[MACX:gnu:*],               [ga_cv_f77_opt="-O3 -funroll-loops"],
[MACX:intel:*],             [ga_cv_f77_opt="-O3 -prefetch -w -cm"],
[NEC64:*:*],                [ga_cv_f77_opt="-Cvsafe -size_t64"],
[NEC:*:*],                  [ga_cv_f77_opt="-Cvsafe"],
[SOLARIS64:fujitsu:*],      [ga_cv_f77_opt="-fw -Kfast -KV9FMADD"],
[SOLARIS64:gnu:*],          [ga_cv_f77_opt="-xs -dalign -xarch=v9"],
[SOLARIS64:gnu:i386],       [ga_cv_f77_opt="-xs -dalign -xarch=amd64"],
[SOLARIS:fujitsu:*],        [ga_cv_f77_opt="-fw -Kfast -KV8PFMADD"],
[SOLARIS:gnu:*],            [ga_cv_f77_opt="-xs -dalign"],
[SOLARIS:gnu:i386],         [ga_cv_f77_opt="-xs -dalign -xarch=sse2"],
                            [ga_cv_f77_opt=])
])])
AC_SUBST([GA_FOPT], [$ga_cv_f77_opt])
])dnl
