AC_DEFUN([ACX_SCALAPACK], [


# Check for scalapack
AC_ARG_WITH([scalapack],
[  --with-scalapack          Specifies the location of the SCALAPACK libraries.], [
case $withval in
  no)
    scalapacklibs=""
    scalapackinc=""
    ;;
  *)
    scalapacklibs="$withval/lib/libscalapack.a"
    ;;
esac
],[ scalapacklibs=""
    scalapackinc=""  ])

CFLAGS="${CXXFLAGS} ${scalapacklibs}"

# Check if we can link against mpi
AC_MSG_CHECKING([if we can link to SCALAPACK])
AC_LANG_PUSH([C])
AC_TRY_LINK( [], [ extern void   Cblacs_pinfo( int* mypnum, int* nprocs); int x, y; Cblacs_pinfo(&x, &y); ],
  [
    HAVE_SCALAPACK=1
    AC_MSG_RESULT([yes])
    AC_DEFINE_UNQUOTED(HAVE_SCALAPACK,${HAVE_SCALAPACK})
    AC_SUBST(HAVE_SCALAPACK)
    AC_DEFINE(HAVE_SCALAPACK)
    AC_SUBST(scalapacklibs)
  ],
  [ HAVE_SCALAPACK=0
    AC_MSG_RESULT([no]) 
  ]
)

AC_LANG_POP([C])

])
