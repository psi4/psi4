AC_DEFUN([ACX_ENABLE_DEBUGGING], [
  acx_enable_debugging="no"
  acx_enable_debugging_flags=""
  AC_ARG_ENABLE([debugging],
    [AC_HELP_STRING([--enable-debugging@<:@=yes|no|OPTION@:>@],
      [Enable debugging C and C++ compilers @<:@default=no@:>@]) ],
    [
      case $enableval in
        yes)
          acx_enable_debugging="yes"
          acx_enable_debugging_flags="-g"
        ;;
        no)
        ;;
        *)
          acx_enable_debugging="yes"
          acx_enable_debugging_flags="-g$enableval"
        ;;
      esac
    ])

  if test $acx_enable_debugging != no; then
    ACX_CHECK_COMPILER_FLAG([C], [CFLAGS], [$acx_enable_debugging_flags],
      [CFLAGS="$CFLAGS $acx_enable_debugging_flags"],
      [AC_MSG_WARN([$CC does not accept $acx_enable_debugging_flags, no debugging flags will be used.])])
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$acx_enable_debugging_flags],
      [CXXFLAGS="$CXXFLAGS $acx_enable_debugging_flags"],
      [AC_MSG_WARN([$CXX does not accept $acx_enable_debugging_flags, no debugging flags will be used.])])
  fi
])