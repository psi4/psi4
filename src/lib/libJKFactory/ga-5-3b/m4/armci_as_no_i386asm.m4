# ARMCI_AS_NO_I386ASM
# -------------------
# On certain systems (Fujitus thus far) it is necessary to set NO_I386ASM to
# avoid compiling inline assembly.
AC_DEFUN([ARMCI_AS_NO_I386ASM], [
AC_CACHE_CHECK([whether NO_I386ASM is needed], [armci_cv_as_no_i386asm],
    [AS_IF([test "x$ax_cv_c_compiler_vendor" = xfujitsu],
        [armci_cv_as_no_i386asm=yes],
        [armci_cv_as_no_i386asm=no])])
AS_IF([test "x$armci_cv_as_no_i386asm" = xyes],
    [AC_DEFINE([NO_I386ASM], [1], [define when inline asm is not supported])])
])dnl
