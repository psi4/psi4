# GA_UNION_SEMUN
# --------------
# Not all packages define union semun, even those that require it.
AC_DEFUN([GA_UNION_SEMUN],
[AC_CACHE_CHECK([for union semun in sys/sem.h], [ga_cv_union_semun],
[AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>]],
        [[union semun arg;
semctl(0, 0, 0, arg);]])],
    [ga_cv_union_semun=yes],
    [ga_cv_union_semun=no])])
AS_IF([test x$ga_cv_union_semun = xyes],
    [AC_DEFINE([HAVE_UNION_SEMUN], [1],
        [define if sys/sem.h has union semun])])
])# GA_UNION_SEMUN
