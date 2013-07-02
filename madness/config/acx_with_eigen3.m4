AC_DEFUN([ACX_WITH_EIGEN3],
[
  acx_with_eigen3=no
  AC_ARG_WITH([eigen3],
    [AS_HELP_STRING([--with-eigen3@<:@=DIR@:>@], [Build with Eigen3 headers.])],
    [
      case $withval in
      yes)
        acx_with_eigen3=yes
      ;;
      no)
        acx_with_eigen3=no
      ;;
      *)
        acx_with_eigen3=yes
        CPPFLAGS="-I$withval $CPPFLAGS"
      ;;
      esac
    ]
  )
  
  if test "$acx_with_eigen3" != no; then
    echo "Eigen3 option required"
    # Check for the pressence of Eigen3 header files.
    AC_CHECK_HEADER([Eigen/Dense], [],
      [AC_MSG_ERROR([Unable to find the Eigen/Dense header file.])])
    AC_DEFINE([MADNESS_HAS_EIGEN3], [1], 
      [Madness will use Eigen3 for linear algebra operations])
  fi
])
