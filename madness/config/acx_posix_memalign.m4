AC_DEFUN([ACX_POSIX_MEMALIGN], [
    AC_CHECK_FUNC([posix_memalign],[gotpm=1], [gotpm=0])
    if test $gotpm = 1; then
        AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [Set if have posix_memalign])
    elif test "$ON_A_MAC" = "yes"; then
        AC_DEFINE([HAVE_POSIX_MEMALIGN], [0], [Set if have posix_memalign])
    else
        AC_MSG_WARN([[   posix_memalign NOT FOUND ... enabling override of new/delete ... THIS WILL BE SLOW ]])
        AC_DEFINE([WORLD_GATHER_MEM_STATS], [1], [Set if MADNESS gathers memory statistics])
    fi

    # look for both version of posix_memalign, with and without throw()
    gotpmh=0
    if test $gotpm = 1; then
        AC_MSG_CHECKING([if missing declaration of posix_memalign in stdlib.h])
        AC_LANG_SAVE
        AC_LANG([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stddef.h>
#include <stdlib.h>
extern "C"  int posix_memalign(void **memptr, size_t alignment, size_t size) throw();]],
[[void *m; posix_memalign(&m, 16, 1024);]])],
         [AC_MSG_RESULT([no])
          gotpmh=1
         ]
        )
        if test $gotpmh = 0; then
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stddef.h>
#include <stdlib.h>
extern "C"  int posix_memalign(void **memptr, size_t alignment, size_t size);]],
[[void *m; posix_memalign(&m, 16, 1024);]])],
           [AC_MSG_RESULT([no])
            gotpmh=1
           ],
           [ AC_DEFINE(MISSING_POSIX_MEMALIGN_PROTO, [1], [Set if the posix_memalign prototype is missing]) 
             AC_MSG_RESULT([yes])
           ]
          )
        fi
        AC_LANG_RESTORE
    fi
])
