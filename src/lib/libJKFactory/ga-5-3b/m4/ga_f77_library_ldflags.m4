# _GA_FC_LIBRARY_LDFLAGS
# ======================
#
# GA HACK!!!
# Until autoconf is fixed to support Sun Studio compilers using -Y, we replace
# _AC_FC_LIBRARY_LDFLAGS.  The rest of the marco is directly copied from
# autoconf-2.67/lib/autoconf/fortran.m4
# ----------------------------------------------------------------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran
# intrinsic and runtime libraries that are required to successfully
# link a Fortran program or shared library.  The output variable
# FLIBS/FCLIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran intrinsic and runtime libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran libraries.  Hence, the macro
# "AC_F77_LIBRARY_LDFLAGS" was created to determine these Fortran
# libraries.
#
# This macro was packaged in its current form by Matthew D. Langston.
# However, nearly all of this macro came from the "OCTAVE_FLIBS" macro
# in "octave-2.0.13/aclocal.m4", and full credit should go to John
# W. Eaton for writing this extremely useful macro.  Thank you John.
AC_DEFUN([_GA_FC_LIBRARY_LDFLAGS],
[_AC_FORTRAN_ASSERT()dnl
_AC_PROG_FC_V
AC_CACHE_CHECK([for _AC_LANG libraries of $[]_AC_FC[]], ac_cv_[]_AC_LANG_ABBREV[]_libs,
[if test "x$[]_AC_LANG_PREFIX[]LIBS" != "x"; then
  ac_cv_[]_AC_LANG_ABBREV[]_libs="$[]_AC_LANG_PREFIX[]LIBS" # Let the user override the test.
else

_AC_PROG_FC_V_OUTPUT

ac_cv_[]_AC_LANG_ABBREV[]_libs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_[]_AC_LANG_ABBREV[]_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
	[[\\/]]*.a | ?:[[\\/]]*.a)
	  _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	      ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
	  ;;
	-bI:*)
	  _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	     [_AC_LINKER_OPTION([$ac_arg], ac_cv_[]_AC_LANG_ABBREV[]_libs)])
	  ;;
	  # Ignore these flags.
	-lang* | -lcrt*.o | -lc | -lgcc* | -lSystem | -libmil | -little \
	  |-LANG:=* | -LIST:* | -LNO:* | -link)
	  ;;
	-lkernel32)
	  test x"$CYGWIN" != xyes && ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
	  ;;
	-[[LRuYz]])
	  # These flags, when seen by themselves, take an argument.
	  # We remove the space between option and argument and re-iterate
	  # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
	     "" | -*);;
	     *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
	  ;;
	-YP,*)
	  for ac_j in `AS_ECHO(["$ac_arg"]) | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
	    _AC_LIST_MEMBER_IF($ac_j, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
			       [ac_arg="$ac_arg $ac_j"
			       ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_j"])
	  done
	  ;;
dnl begin GA addition
	-Y*)
	  for ac_j in `AS_ECHO(["$ac_arg"]) | sed -e 's/-Y/-L/;s/"//g;s/:/ -L/g'`; do
	    _AC_LIST_MEMBER_IF($ac_j, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
			       [ac_arg="$ac_arg $ac_j"
			       ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_j"])
	  done
	  ;;
dnl end GA addition
	-[[lLR]]*)
	  _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
			     ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
	  ;;
	-zallextract*| -zdefaultextract)
	  ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
	  ;;
	  # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`AS_ECHO(["$ac_[]_AC_LANG_ABBREV[]_v_output"]) |
			sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
	_AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_[]_AC_LANG_ABBREV[]_libs)
      ;;
esac
fi # test "x$[]_AC_LANG_PREFIX[]LIBS" = "x"
])
[]_AC_LANG_PREFIX[]LIBS="$ac_cv_[]_AC_LANG_ABBREV[]_libs"
AC_SUBST([]_AC_LANG_PREFIX[]LIBS)
])# _GA_FC_LIBRARY_LDFLAGS


# GA_F77_LIBRARY_LDFLAGS
# ----------------------
# Wrap AC_F77_LIBRARY_LDFLAGS in case user disables Fortran 77.
# When mixing gcc and ifort, we sometimes need to add -lgcc_s to the FLIBS.
# Lastly, Sun Studio compilers are not (yet) completely supported.
AC_DEFUN([GA_F77_LIBRARY_LDFLAGS], [
    GA_MPI_UNWRAP_PUSH()
    AC_F77_LIBRARY_LDFLAGS
    GA_MPI_UNWRAP_POP()
    AC_CACHE_CHECK([whether FLIBS needs -lgcc_s], [ga_cv_flibs_gcc_s],
        [happy=yes
         ga_save_LIBS="$LIBS";  LIBS="$LIBS $FLIBS"
         ga_save_FLIBS="$FLIBS"
         AC_LANG_PUSH([C])
         AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([],[])],
            [ga_cv_flibs_gcc_s=no],
            [LIBS="$LIBS -lgcc_s"
             AC_LINK_IFELSE(
                [AC_LANG_PROGRAM([],[])],
                [FLIBS="$FLIBS -lgcc_s"
                 ga_cv_flibs_gcc_s=yes],
                [happy=no
                 FLIBS="$ga_save_FLIBS"
                 ga_cv_flibs_gcc_s=no])])
         LIBS="$ga_save_LIBS"
         AC_LANG_POP([C])])
dnl hack to support Sun Studio -Y linker flag
    AS_IF([test "x$happy" = xno],
        [AS_UNSET([ac_cv_f77_libs])
         AS_UNSET([FLIBS])
         _GA_FC_LIBRARY_LDFLAGS
         ga_save_LIBS="$LIBS";  LIBS="$LIBS $FLIBS"
         AC_LANG_PUSH([C])
         AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([],[])],
            [happy=yes], [happy=no])
         LIBS="$ga_save_LIBS"
         AC_LANG_POP([C])])
    AS_IF([test "x$happy" = xno],
        [AC_MSG_WARN([FLIBS does not work])])
])dnl
