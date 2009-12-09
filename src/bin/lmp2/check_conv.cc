/*! \file
    \ingroup LMP2
    \brief check the LMP2 convergence
*/
#include "mpi.h"
#include <iostream>
#include <fstream>              // file I/O support
#include <cstdlib>              // support for exit()
#include <cstring>
#include <memory.h>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <psitypes.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;

void LMP2::check_conv() {

  int i, j, a, b, v, ij;
  int *ij_owner, *ij_local, count;
  double rms, rms_sum;

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();


  //   ****  Compute the Energy Difference ****
  DEmp2 = Emp2 - Emp2_old;

  if (iter <= it_diis || diis == 0) {
    //   ****  Compute the RMS for the New and Old Amplitudes  ****
    rms = 0.0;
    v=0;
    count = 0;
    for(i=0, ij=0; i < nocc; i++) {
      for(j=0; j <= i; j++, ij++, v++) {
        if(v%nprocs != myid)
          continue;

        for(a=0; a < pairdom_len[ij]; a++) {
          for(b=0; b < pairdom_len[ij]; b++) {
            if(fabs(T[div][ij_local[ij]][a][b]) > 1e-14 && fabs(T[dmat1][ij_local[ij]][a][b]) > 1e-14) {
              if(i != j)
                rms += 2 * (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]) *
                       (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]);
              else
                rms += (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]) *
                       (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]);
            }
            count++;
          }
        }
      }
    }
  }
  else {
    //   ****  Compute the RMS for the New and Old Amplitudes  ****
    rms = 0.0;
    v=0;
    count = 0;
    for(i=0, ij=0; i < nocc; i++) {
      for(j=0; j <= i; j++, ij++, v++) {
        if(v%nprocs != myid)
          continue;

        for(a=0; a < pairdom_len[ij]; a++) {
          for(b=0; b < pairdom_len[ij]; b++) {
            if(fabs(T[div][ij_local[ij]][a][b]) > 1e-14 && fabs(T[dmat1][ij_local[ij]][a][b]) > 1e-14) {
              if(i != j)
                rms += 2 * (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]) *
                       (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]);
              else
                rms += (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]) *
                       (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]);
            }
            count++;
          }
        }
      }
    }
  }

  Drms = 0.0;
  MPI_Allreduce(&rms, &Drms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Drms = sqrt(Drms/count);

  if(fabs(DEmp2) < econv && fabs(Drms) < rmsconv || iter >= maxiter) {
    conv = 1;
    if(iter >= maxiter)
      if(myid == 0)
        fprintf(outfile, "LMP2 has not converged in the maximum number of iterations.\n maxiter = %d\n", maxiter);
  }
  else conv = 0;

}

}} // namespace psi::lmp2
