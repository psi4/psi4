AC_DEFUN([ACX_WITH_LIBUNWIND], [
  acx_with_libunwind=""
  AC_ARG_WITH([libunwind],
    [AS_HELP_STRING([--with-libunwind@<:@=Install DIR@:>@],
      [Enables use of libunwind, required for Google profiler and heap checker])],
    [
      case $withval in
      yes)
        acx_with_libunwind="yes"
      ;;
      no)
        acx_with_libunwind="no"
      ;;
      *)
        CPPFLAGS="$CPPFLAGS -I$with_libunwind/include"
        LIBS="$LIBS -L$with_libunwind/lib"
        acx_with_libunwind="$withval"
      esac
    ],
    [acx_with_libunwind="no"]
  )
  if test $acx_with_libunwind != "no"; then
    AC_LANG_SAVE
    AC_LANG_C
    AC_CHECK_HEADER([libunwind.h], [], [AC_MSG_ERROR([Unable to compile with the libunwind.])])
    AC_CHECK_LIB([unwind], [_Unwind_GetRegionStart], [LIBS="$LIBS -lunwind"], [AC_MSG_ERROR(["Unable to link with libunwind])])
    AC_DEFINE(MADNESS_HAS_LIBUNWIND, [1], [Define if should use libunwind for Google performance tools])
    AC_LANG_RESTORE
  fi
])
