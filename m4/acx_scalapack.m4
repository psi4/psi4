AC_DEFUN([ACX_SCALAPACK], [


# Check for scalapack
AC_ARG_WITH([scalapack],
[  --with-scalapack          Specifies the location of the SCALAPACK libraries.], [
case $withval in
  no)
    scalpacklibs=""
    scalpackinc=""
    ;;
  *)
    scalpacklibs="$withval/lib/libscalapack.a"
    ;;
esac
],[ scalpacklibs=""
    scalpackinc=""  ])

CFLAGS="${CXXFLAGS} ${scalpacklibs}"

# Check if we can link against mpi
AC_MSG_CHECKING([if we can link to SCALAPACK])
AC_LANG(C)
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

])
