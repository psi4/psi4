# _GA_F2C_CMDARGS_CONFTEST
# ------------------------
# Link a Fortran 77 program that uses the intrinsics GETARG and IARGC.
AC_DEFUN([_GA_F2C_CMDARGS_CONFTEST],
[AC_LINK_IFELSE([[      program main
      $ga_fxx_module
      integer i, j, l, ier
      character*20 s
      $ga_f77_getarg_decl
      i = 0
      call $ga_f77_getarg($ga_f77_getarg_args)
      i=$ga_f77_iargc()
      if (i .gt. 1) then
          j = i - $ga_f77_iargc()
          j = 1.0 / j
      endif
      end]],
    [],
    [ga_f77_getarg=])
])dnl


# _GA_F2C_CMDARGS
# --------------
# Determine how to access the Fortran 77 command line from C.
#
# Output Effects:
#  The following variables are set:
#    ga_f77_getarg      - Statement to get an argument i into string s
#    ga_f77_getarg_args - Arguments to getarg
#    ga_f77_iargc       - Routine to return the number of arguments
#    ga_fxx_module      - Module command when using Fortran 90 compiler
#    ga_f77_getarg_decl - Declaration of routine used for ga_f77_getarg
# If 'ga_f77_getarg' has a value, then that value and the values for these
# other symbols will be tested first. If no approach is found, all of these
# variables will have empty values.
# 'AC_SUBST' is called for all seven variables.
AC_DEFUN([_GA_F2C_CMDARGS],
[AC_LANG_PUSH([Fortran 77])
dnl AC_ARG_VAR([F77_GETARG_DECL], [Declaration of routine e.g. external GETARG])
dnl AC_ARG_VAR([F77_GETARG],      [Name of routine e.g. getarg, pxfgetarg])
dnl AC_ARG_VAR([F77_GETARG_ARGS], [Arguments to getarg e.g. i,s or i,s,l,ier])
dnl AC_ARG_VAR([F77_IARGC],       [Name of routine e.g. iargc, ipxfargc])
dnl AC_ARG_VAR([FXX_MODULE],
dnl     [Module command when using Fortran 90 compiler e.g. use f90_unix])
AC_MSG_CHECKING([for routines to access the command line from Fortran])
# User-specified values
dnl AS_IF([test "x$ga_f77_getarg" = x],
dnl     [ga_fxx_module=$FXX_MODULE
dnl     ga_f77_getarg_decl=$F77_GETARG_DECL
dnl     ga_f77_getarg=$F77_GETARG
dnl     ga_f77_getarg_args=$F77_GETARG_ARGS
dnl     ga_f77_iargc=$F77_IARGC
dnl     _GA_F2C_CMDARGS_CONFTEST()])
# Standard practice, uppercase
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module=
    ga_f77_getarg_decl="external GETARG"
    ga_f77_getarg=GETARG
    ga_f77_getarg_args="i,s"
    ga_f77_iargc=IARGC
    _GA_F2C_CMDARGS_CONFTEST()])
# Standard practice, lowercase
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module=
    ga_f77_getarg_decl="external getarg"
    ga_f77_getarg=getarg
    ga_f77_getarg_args="i,s"
    ga_f77_iargc=iargc
    _GA_F2C_CMDARGS_CONFTEST()])
# Posix alternative
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module=
    ga_f77_getarg_decl="external pxfgetarg"
    ga_f77_getarg=pxfgetarg
    ga_f77_getarg_args="i,s,l,ier"
    ga_f77_iargc=ipxfargc
    _GA_F2C_CMDARGS_CONFTEST()])
# Nag f90_unix_env module
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module="        use f90_unix_env"
    ga_f77_getarg_decl=
    ga_f77_getarg=getarg
    ga_f77_getarg_args="i,s"
    ga_f77_iargc=iargc
    _GA_F2C_CMDARGS_CONFTEST()])
# Nag f90_unix module
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module="        use f90_unix"
    ga_f77_getarg_decl=
    ga_f77_getarg=getarg
    ga_f77_getarg_args="i,s"
    ga_f77_iargc=iargc
    _GA_F2C_CMDARGS_CONFTEST()])
# gfortran won't find getarg if it is marked as external
AS_IF([test "x$ga_f77_getarg" = x],
    [ga_fxx_module=
    ga_f77_getarg_decl="intrinsic GETARG"
    ga_f77_getarg=GETARG
    ga_f77_getarg_args="i,s"
    ga_f77_iargc=IARGC
    _GA_F2C_CMDARGS_CONFTEST()])
AS_IF([test "x$ga_f77_getarg" = x],
    [AC_MSG_RESULT([no])
    AC_MSG_ERROR([Could not find way to access Fortran cmd line from C])])
AC_MSG_RESULT([yes])
AC_DEFINE_UNQUOTED([F77_GETARG_DECL], [$ga_f77_getarg_decl],
    [Declaration of routine e.g. external GETARG])
AC_DEFINE_UNQUOTED([F77_GETARG], [$ga_f77_getarg],
    [Name of routine e.g. getarg, pxfgetarg])
AC_DEFINE_UNQUOTED([F77_GETARG_ARGS], [$ga_f77_getarg_args],
    [Arguments to getarg e.g. i,s or i,s,l,ier])
AC_DEFINE_UNQUOTED([F77_IARGC], [$ga_f77_iargc],
    [Name of routine e.g. iargc, ipxfargc])
AC_DEFINE_UNQUOTED([FXX_MODULE], [$ga_fxx_module],
    [Module command when using Fortran 90 compiler e.g. use f90_unix])
AC_F77_FUNC([F2C_GETARG])
AC_F77_FUNC([F2C_IARGC])
AC_LANG_POP([Fortran 77])
])dnl

# GA_F2C_CMDARGS
# --------------
# Determine how to access the Fortran 77 command line from C.
AC_DEFUN([GA_F2C_CMDARGS], [
_GA_F2C_CMDARGS
AC_SUBST([F2C_GETARG])
AC_SUBST([F2C_IARGC])
])dnl
