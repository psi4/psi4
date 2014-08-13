# TCGMSG_ENABLE_TIMINGS
# ---------------------
# This was only defined for LINUX[64] previously. Does it do any harm to
# define it permanently?
AC_DEFUN([TCGMSG_ENABLE_TIMINGS],
[AC_DEFINE([TCGMSG_TIMINGS], [1], [Gather timing information for TCGMSG])
])dnl
