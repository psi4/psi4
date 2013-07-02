AC_DEFUN([ACX_MPI], [
  HAVE_MPI="0"
  AC_CXX_PROCESS_CHECK([mpi.h],[#include <mpi.h>],,[HAVE_MPI="1"],[HAVE_MPI="0"])
  AC_MSG_CHECKING([if we can link to MPI])
  AC_TRY_LINK([#include <mpi.h>],[ MPI_Finalize(); ],
    [HAVE_MPI="1"; AC_MSG_RESULT(yes)],
    [HAVE_MPI="0"; AC_MSG_RESULT(no)])

  if test $HAVE_MPI == "1"; then
    AC_DEFINE_UNQUOTED(HAVE_MPI,1)
    AC_SUBST(HAVE_MPI)
  fi
])
