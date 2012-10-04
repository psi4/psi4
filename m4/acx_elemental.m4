AC_DEFUN([ACX_ELEMENTAL], [


# Check for elemental
AC_ARG_WITH([elemental],
[  --with-elemental          Specifies the location of the ELEMENTAL libraries.], [
case $withval in
  no)
    elemlibs=""
    eleminc=""
    ;;
  *)
    elemlibs="$withval/lib/libelemental.a $withval/lib/libplcg.a $withval/lib/libpmrrr.a $withval/lib/liblapack-addons.a"
    eleminc="-I$withval/include"
    ;;
esac
],[ elemlibs=""
    eleminc=""  ])

CXXFLAGS="${CXXFLAGS} ${elemlibs} ${eleminc}"

# Check if we can link against mpi
AC_MSG_CHECKING([if we can link to ELEMENTAL])
AC_LANG_PUSH([C++])
AC_TRY_LINK( [#include <elemental.hpp>], [ elem::Finalize(); ],
  [
    HAVE_ELEMENTAL=1
    AC_MSG_RESULT([yes])
    AC_DEFINE_UNQUOTED(HAVE_ELEMENTAL,${HAVE_ELEMENTAL})
    AC_SUBST(HAVE_ELEMENTAL)
    AC_DEFINE(HAVE_ELEMENTAL)
    AC_SUBST(elemlibs)
  ],
  [ HAVE_ELEMENTAL=0
    USE_PSIDISTMATRIX=1
    AC_DEFINE_UNQUOTED(USE_PSIDISTMATRIX,${USE_PSIDISTMATRIX})
    AC_SUBST(USE_PSIDISTMATRIX)
    AC_MSG_RESULT([no]) 
  ]
)

AC_LANG_POP([C++])

])
