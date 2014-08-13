# COMEX_ARG_PARSE(ARG, VAR_LIBS, VAR_LDFLAGS, VAR_CPPFLAGS)
# ------------------------------------------------------
# Parse whitespace-separated ARG into appropriate LIBS, LDFLAGS, and
# CPPFLAGS variables.
# some examples (not exhaustive):
# $arg="/apps/mpi/hp/2.03.01.00/include/64 /apps/mpi/hp/2.03.01.00/lib/linux_amd64 -lmtmpi"
# $arg="/usr/local"
# $arg="-I/usr/local/include -L/usr/local/lib -lsomelib"
# $arg="/usr/lib/mpich-shmem -lfmpich-shmem -lmpich-shmem -lpmpich-shmem"
# $arg="/usr/lib/mpich-shmem/include /usr/lib/mpich-shmem/lib -lfmpich-shmem -lmpich-shmem -lpmpich-shmem"
# $arg="/usr/local/lib64 /usr/local/include64"
#
# Special Cases:
# -mkl[=arg]    Intel compiler with special -mkl flag for headers and linking
AC_DEFUN([COMEX_ARG_PARSE],
[AC_COMPUTE_INT([comex_arg_parse_sizeof_voidp], [(long int) (sizeof (void*))])
for arg in $$1 ; do
    AS_CASE([$arg],
        [yes],          [],
        [no],           [],
        [-l*],          [$2="$$2 $arg"],
        [-L*],          [$3="$$3 $arg"],
        [-WL*],         [$3="$$3 $arg"],
        [-Wl*],         [$3="$$3 $arg"],
        [-I*],          [$4="$$4 $arg"],
        [*.a],          [$2="$$2 $arg"],
        [*.so],         [$2="$$2 $arg"],
        [*lib],         [AS_IF([test -d $arg], [$3="$$3 -L$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*lib/],        [AS_IF([test -d $arg], [$3="$$3 -L$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*lib64],       [AS_IF([test -d $arg], [$3="$$3 -L$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*lib64/],      [AS_IF([test -d $arg], [$3="$$3 -L$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*include],     [AS_IF([test -d $arg], [$4="$$4 -I$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*include/],    [AS_IF([test -d $arg], [$4="$$4 -I$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*include64],   [AS_IF([test -d $arg], [$4="$$4 -I$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [*include64/],  [AS_IF([test -d $arg], [$4="$$4 -I$arg"],
                            [AC_MSG_WARN([$arg of $1 not parsed])])],
        [-mkl*],        [$2="$$2 $arg"],
        [comex_arg_parse_ok=no])
# $arg didn't fit the most common cases
# check for subdirectories e.g. lib,include
    AS_IF([test "x$comex_arg_parse_ok" = xno],
        [AS_IF([test "x$comex_arg_parse_sizeof_voidp" = x8],
            [AS_IF([test -d $arg/lib64],    [$3="$$3 -L$arg/lib64"; comex_arg_parse_ok=yes],
                   [test -d $arg/lib],      [$3="$$3 -L$arg/lib"; comex_arg_parse_ok=yes])
             AS_IF([test -d $arg/include64],[$4="$$4 -I$arg/include64"; comex_arg_parse_ok=yes],
                   [test -d $arg/include],  [$4="$$4 -I$arg/include"; comex_arg_parse_ok=yes])],
            [AS_IF([test -d $arg/lib],      [$3="$$3 -L$arg/lib"; comex_arg_parse_ok=yes])
             AS_IF([test -d $arg/include],  [$4="$$4 -I$arg/include"; comex_arg_parse_ok=yes])])])
# $arg still unknown, look for "lib" and "include" anywhere...
    AS_IF([test "x$comex_arg_parse_ok" = xno],
        [AS_CASE([$arg],
            [*lib*],    [AS_IF([test -d $arg], [$3="$$3 -L$arg"; comex_arg_parse_ok=yes])],
            [*include*],[AS_IF([test -d $arg], [$4="$$4 -I$arg"; comex_arg_parse_ok=yes])])])
# warn user that $arg fell through
     AS_IF([test "x$comex_arg_parse_ok" = xno],
        [AC_MSG_WARN([$arg of $1 not parsed])])
done])dnl
