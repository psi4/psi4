# COMEX_WITH_MPI
# --------------
# Establishes all things related to the Message Passing Interface.
# This includes the compilers to use (either standard or MPI wrappers)
# or the proper linker flags (-L), libs (-l) or preprocessor directives (-I).
# Yes, it's a beefy AC macro, but because when MPI is desired it replaces the
# usual compiler the order here is necessary and it is all interdependent.
AC_DEFUN([COMEX_WITH_MPI], [
# MPI_* vars might exist in environment, but they are really internal.
# Reset them.
MPI_LIBS=
MPI_LDFLAGS=
MPI_CPPFLAGS=
AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi[[=ARG]]],
        [path to MPI; leave ARG blank to use MPI compiler wrappers in PATH])],
    [],
    [with_mpi=yes])
AS_IF([test "x$with_mpi" = xyes],
    [with_mpi_wrappers=yes],
    [with_mpi_wrappers=no])
dnl postpone parsing with_mpi until we know sizeof(void*)
dnl AS_IF([test x$with_mpi_wrappers = xno],
dnl     [COMEX_ARG_PARSE([with_mpi], [MPI_LIBS], [MPI_LDFLAGS], [MPI_CPPFLAGS])])
AC_SUBST([MPI_LIBS])
AC_SUBST([MPI_LDFLAGS])
AC_SUBST([MPI_CPPFLAGS])
])dnl
