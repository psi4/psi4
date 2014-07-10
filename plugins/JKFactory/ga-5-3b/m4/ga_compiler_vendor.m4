# GA_COMPILER_VENDOR
# ------------------
# other compilers may define __GNUC__, like intel or pathscale
# that's why we put it later in the for loop, like we do with our compiler
# checks since GCC is so pervasive
AC_DEFUN([GA_COMPILER_VENDOR], [
AS_VAR_PUSHDEF([ga_cv_compiler_vendor],
               [ga_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])
AC_CACHE_CHECK([for _AC_LANG compiler vendor], [ga_cv_compiler_vendor], [
ga_save_ac_ext="$ac_ext"
AC_LANG_CASE([Fortran],    [ac_ext=F])
AC_LANG_CASE([Fortran 77], [ac_ext=F])
ga_cv_compiler_vendor=unknown
ga_cpp_vendor_symbols=
for vendor in intel ibm pathscale amd cray gnu sun hp dec borland comeau kai lcc metrowerks sgi microsoft watcom portland fujitsu
do
AS_CASE([$vendor],
[amd],       [ga_cpp_vendor_symbols="defined(__OPEN64__)"],
[borland],   [ga_cpp_vendor_symbols="defined(__BORLANDC__) || defined(__TURBOC__)"],
[comeau],    [ga_cpp_vendor_symbols="defined(__COMO__)"],
[cray],      [ga_cpp_vendor_symbols="defined(_CRAYC) || defined(_ADDR64)"],
[dec],       [ga_cpp_vendor_symbols="defined(__DECC) || defined(__DECCXX) || defined(__DECC_VER) || defined(__DECCXX_VER)"],
[fujitsu],   [ga_cpp_vendor_symbols="defined(__fcc__) || defined(__fcc_version__) || defined(_FCC_VER) || defined(__FCC_VER_)"],
[gnu],       [ga_cpp_vendor_symbols="defined(__GNUC__)"],
[hp],        [ga_cpp_vendor_symbols="defined(__HP_cc) || defined(__HP_aCC)"],
[ibm],       [ga_cpp_vendor_symbols="defined(__xlc__) || defined(__xlC__) || defined(__IBMC__) || defined(__IBMCPP__)"],
[intel],     [ga_cpp_vendor_symbols="defined(__ICC) || defined(__ECC) || defined(__INTEL_COMPILER)"],
[kai],       [ga_cpp_vendor_symbols="defined(__KCC)"],
[lcc],       [ga_cpp_vendor_symbols="defined(__LCC__)"],
[metrowerks],[ga_cpp_vendor_symbols="defined(__MWERKS__)"],
[microsoft], [ga_cpp_vendor_symbols="defined(_MSC_VER)"],
[pathscale], [ga_cpp_vendor_symbols="defined(__PATHCC__) || defined(__PATHSCALE__)"],
[portland],  [ga_cpp_vendor_symbols="defined(__PGI)"],
[sgi],       [ga_cpp_vendor_symbols="defined(__sgi) || defined(sgi)"],
[sun],       [ga_cpp_vendor_symbols="defined(__SUNPRO_C) || defined(__SUNPRO_CC)"],
[watcom],    [ga_cpp_vendor_symbols="defined(__WATCOMC__)"])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#if !($ga_cpp_vendor_symbols)
    chokeonthis
#endif
])], [ga_cv_compiler_vendor=$vendor; break])
done
ga_cpp_vendor_symbols=
ac_ext="$ga_save_ac_ext"
])
AS_VAR_POPDEF([ga_cv_compiler_vendor])
])dnl
