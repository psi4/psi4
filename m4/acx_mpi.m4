AC_DEFUN([ACX_MPI], [
  HAVE_MPI="0"
  # Check for an MPI C compiler
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)

  # Check for an MPI C++ compiler
  AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
  AC_CHECK_PROGS(MPICXX, mpicxx mpic++ mpiCC mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
 
  # set have_mpi 
  if test $MPICC != $CC; then
    if test $MPICXX != $CXX; then
      HAVE_MPI="1"
      AC_DEFINE_UNQUOTED(HAVE_MPI,1)
    else
      AC_MSG_WARN([MPICXX was not found.])
    fi
  else
    AC_MSG_WARN([MPICC was not found.])
  fi
])
