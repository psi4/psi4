AC_DEFUN([ACX_ENABLE_OPTIMIZATION], [
  default_optimization=""
  case $CXXVENDOR in
    GNU)
      default_optimization="3"
    ;;
    clang)
      default_optimization="3"
    ;;
    Pathscale)
      default_optimization="fast"
    ;;
    Portland)
      default_optimization="3"
    ;;
    Intel)
      default_optimization="3"
    ;;
    IBM)
      default_optimization="3"
    ;;
    *)
      default_optimization="2"
    ;;
  esac
  acx_enable_optimization="yes"
  acx_enable_optimization_flags=""
  AC_ARG_ENABLE([optimization],
    [AC_HELP_STRING([--enable-optimization@<:@=yes|no|OPTION@:>@],
      [Enable optimization for C and C++ @<:@default=yes@:>@]) ],
    [
      case $enableval in
        yes)
          acx_enable_optimization_flags="-O$default_optimization"
        ;;
        no)
          acx_enable_optimization="no"
        ;;
        *)
          acx_enable_optimization_flags="-O$enableval"
        ;;
      esac
    ],
    [
      if test "x$acx_enable_debugging" == xno; then
        acx_enable_optimization_flags="-O$default_optimization"
      else
        acx_enable_optimization_flags="-O0"
        AC_MSG_WARN([Optimizations is disabled, because debugging is enabled. Add --enable-optimization to overide this behavior.])
      fi
    ]
  )
  
  if test $acx_enable_optimization != no; then
    ACX_CHECK_COMPILER_FLAG([C], [CFLAGS], [$acx_enable_optimization_flags],
      [CFLAGS="$CFLAGS $acx_enable_optimization_flags"],
      [AC_MSG_WARN([$CC does not accept $acx_enable_optimization_flags, no optimization flags will be used.])])
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$acx_enable_optimization_flags],
      [CXXFLAGS="$CXXFLAGS $acx_enable_optimization_flags"],
      [AC_MSG_WARN([$CXX does not accept $acx_enable_optimization_flags, no optimization flags will be used.])])
  fi

])