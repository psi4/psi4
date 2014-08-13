# GA_CXX_OPT()
# ------------
# Determine TARGET-/compiler-specific CXXFLAGS for optimization.
AC_DEFUN([GA_CXX_OPT], [
AC_REQUIRE([GA_TARGET64])
AC_REQUIRE([GA_ENABLE_OPT])
AC_ARG_VAR([GA_CXXOPT], [GA C++ optimization flags])
AC_CACHE_CHECK([for specific C++ optimizations], [ga_cv_cxx_opt], [
AS_IF([test "x$GA_CXXOPT" != x], [ga_cv_cxx_opt="$GA_CXXOPT"], [ga_cv_cxx_opt=])
AS_IF([test "x$ga_cv_cxx_opt" = x && test "x$enable_opt" = xyes], [
AS_CASE([$ga_cv_target:$ga_cv_cxx_compiler_vendor:$host_cpu],
[LINUX:*:*],                [ga_cv_cxx_opt="-O0"],
[NEC64:*:*],                [ga_cv_cxx_opt="-Cvsafe -size_t64"],
[NEC:*:*],                  [ga_cv_cxx_opt="-Cvsafe"],
                            [ga_cv_cxx_opt=])
])])
AC_SUBST([GA_CXXOPT], [$ga_cv_cxx_opt])
])dnl
