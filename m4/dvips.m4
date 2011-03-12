dnl
dnl Check for dvips
dnl
AC_DEFUN([AC_PROG_DVIPS],
[AC_PROVIDE([$0])
AC_CHECK_PROG(DVIPS, dvips, dvips)]
)dnl

dnl
dnl Check for ps2pdf
dnl
AC_DEFUN([AC_PROG_PS2PDF],
[AC_PROVIDE([$0])
AC_CHECK_PROG(PS2PDF, ps2pdf, ps2pdf)]
)dnl
