# ARMCI_CXX_OPT()
# ---------------
# Determine TARGET-/compiler-specific CXXFLAGS and FFLAGS for optimization.
AC_DEFUN([ARMCI_CXX_OPT], [
AC_REQUIRE([GA_TARGET64])
AC_REQUIRE([GA_ENABLE_OPT])
AC_REQUIRE([GA_ARMCI_NETWORK])
AC_ARG_VAR([ARMCI_CXXOPT], [ARMCI C++ optimization flags])
AC_CACHE_CHECK([for specific C++ optimizations], [armci_cv_cxx_opt], [
AS_IF([test "x$ARMCI_CXXOPT" != x], [armci_cv_cxx_opt="$ARMCI_CXXOPT"], [armci_cv_cxx_opt=])
AS_IF([test "x$armci_cv_cxx_opt" = x && test "x$enable_opt" = xyes], [
AS_CASE([$ga_cv_target:$ga_cv_cxx_compiler_vendor:$host_cpu:$ga_armci_network],
[LINUX:*:*:*],              [armci_cv_cxx_opt="-O0"],
[NEC64:*:*:*],              [armci_cv_cxx_opt="-Cvsafe -size_t64"],
[NEC:*:*:*],                [armci_cv_cxx_opt="-Cvsafe"],
                            [armci_cv_cxx_opt=])
])])
AC_SUBST([ARMCI_CXXOPT],    [$armci_cv_cxx_opt])
])dnl
