AC_DEFUN([ACX_WITH_TBB], [
  acx_with_tbb_include="no"
  acx_with_tbb_lib="no"
  acx_with_tbb="no"
  
  # Configure madness to use Intel TBB and specify the include path.
  AC_ARG_WITH([tbb-include],
    [AS_HELP_STRING([--with-tbb-include@<:@=DIR@:>@],
      [Enables use of Intel TBB as the task scheduler.])],
    [
      case $withval in
      yes)
        AC_MSG_ERROR([You must specify a directory for --with-tbb-include.])
      ;;
      no)
      ;;
      *)
        CPPFLAGS="$CPPFLAGS -I$withval"
        acx_with_tbb_include="yes"
        acx_with_tbb="yes"
      esac
    ]
  )
  
  
  # Configure madness to use Intel TBB and specify the library path.
  AC_ARG_WITH([tbb-lib],
    [AS_HELP_STRING([--with-tbb-lib@<:@=DIR@:>@],
      [Enables use of Intel TBB as the task scheduler.])],
    [
      case $withval in
      yes)
        AC_MSG_ERROR([You must specify a directory for --with-tbb-lib.])
      ;;
      no)
      ;;
      *)
        LIBS="$LIBS -L$withval"
        acx_with_tbb_lib="yes"
        acx_with_tbb="yes"
      esac
    ]
  )
  
  # Configure madness to use Intel TBB
  AC_ARG_WITH([tbb],
    [AS_HELP_STRING([--with-tbb@<:@=Install DIR@:>@],
      [Enables use of Intel TBB as the task scheduler.])],
    [
      case $withval in
      yes)
        acx_with_tbb="yes"
      ;;
      no)
      ;;
      *)
        if test "$acx_with_tbb_include" == no; then
          CPPFLAGS="$CPPFLAGS -I$withval/include"
        fi
        if test "$acx_with_tbb_lib" == no; then
          LIBS="$LIBS -L$withval/lib"
        fi
        acx_with_tbb="yes"
      esac
    ],
    [acx_with_tbb="no"]
  )
  
  # Check that we can compile with Intel TBB
  if test $acx_with_tbb != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    
    # Check for Intel TBB header.
    AC_CHECK_HEADER([tbb/tbb.h], [], [AC_MSG_ERROR([Unable to compile with Intel TBB.])])
    
    # Check for Intel TBB library.
    if test "x$acx_enable_debugging" == xno; then
      AC_CHECK_LIB([tbb], [TBB_runtime_interface_version], [LIBS="$LIBS -ltbb"], [AC_MSG_ERROR(["Unable to link with Intel TBB])])
    else
      AC_CHECK_LIB([tbb_debug], [TBB_runtime_interface_version], [LIBS="$LIBS -ltbb_debug"], [AC_MSG_ERROR([Unable to link with Intel TBB.])])
      CPPFLAGS="$CPPFLAGS -DTBB_USE_DEBUG=1"
      AC_MSG_WARN([Linking with the debug variant of Intel TBB.])
    fi
    
    AC_DEFINE(MADNESS_HAS_TBB, [1], [Define if Intel TBB is available.])
    AC_LANG_RESTORE
  fi
])