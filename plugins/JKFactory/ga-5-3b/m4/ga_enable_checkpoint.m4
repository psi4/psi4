# GA_ENABLE_CHECKPOINT
# --------------------
# Whether to enable checkpoinging. AC_DEFINEs ENABLE_CHECKPOINT.
AC_DEFUN([GA_ENABLE_CHECKPOINT],
[AC_ARG_ENABLE([checkpoint],
    [AS_HELP_STRING([--enable-checkpoint], [enable checkpointing])],
    [enable_checkpoint=yes
    AC_DEFINE([ENABLE_CHECKPOINT], [1], [Define if checkpointing is enabled])],
    [enable_checkpoint=no])
AM_CONDITIONAL([ENABLE_CHECKPOINT], [test x$enable_checkpoint = xyes])
])dnl
