AC_DEFUN([ACX_MPI], [
        # We were using the full macro from ACX but this forced AC_PROG_C or AC_PROG_CXX to 
        # be run before we had overridden the compilers which meant that some confdef.h
        # entries were incorrect (specifically std::exit problem with PGI)

        AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
        acx_mpi_save_CC="$CC"
        CC="$MPICC"
        AC_SUBST(MPICC)

        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        AC_CHECK_PROGS(MPICXX, mpicxx mpic++ mpiCC mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)

        #AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
        #AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, $F77)
        #acx_mpi_save_F77="$F77"
        #F77="$MPIF77"
        #AC_SUBST(MPIF77)
])
