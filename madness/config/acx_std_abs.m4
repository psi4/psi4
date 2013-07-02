AC_DEFUN([ACX_STD_ABS], [
    AC_MSG_CHECKING([std::abs(long)])
    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <cmath>
#include <cstdlib>
long (*absptr)(long) = &std::abs; 
long a = -1;  
long b = std::abs(a);
]])],
                   [AC_DEFINE(HAVE_STD_ABS_LONG,[1],[Define if have std::abs(long)]) have_abs_long=yes],
                   [have_abs_long=no])
    AC_MSG_RESULT([$have_abs_long])
    if test X"$have_abs_long" = Xno; then
        AC_MSG_CHECKING([std::labs(long)])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <cmath>
#include <cstdlib>
long (*labsptr)(long) = &std::labs; 
long a = -1l;  
long b = labs(a);
]])],
                       [AC_DEFINE(HAVE_LABS,[1],[Define if have std::labs(long)]) have_std_labs=yes],
                       [have_std_labs=no])
        AC_MSG_RESULT([$have_std_labs])
    fi
])
