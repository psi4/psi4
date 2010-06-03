/*! \file
    \ingroup LMP2
    \brief check the LMP2 convergence
*/
//#include "mpi.h"
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdlib>              // support for exit()
//#include <cstring>
//#include <memory.h>
//#include <assert.h>
#include <math.h>
//#include <cmath>
//#include <sstream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <psitypes.h>
//#include <libipv1/ip_lib.h>
//#include <libpsio/psio.h>
//#include <libciomr/libciomr.h>
//#include <libiwl/iwl.h>
//#include <psifiles.h>
//#include <libchkpt/chkpt.h>
//#include <libqt/qt.h>
//#include <libchkpt/chkpt.hpp>
//#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>
#include <libmints/mints.h>
#define EXTERN
#include "globals.h"

namespace psi{

namespace lmp2{

void LMP2::check_conv() {

  int i, j, a, b, ij;
  int *ij_owner, *ij_local, count;
  int **ij_map;
  int **pairdomain, *pairdom_len;

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();
  ij_map = get_ij_map();
  pairdomain = compute_pairdomain(ij_map);
  pairdom_len = compute_pairdomlen(ij_map);

  //   ****  Compute the Energy Difference ****
  DEmp2 = Emp2 - Emp2_old;

  if (iter <= it_diis || diis == 0) {
    //   ****  Compute the RMS for the New and Old Amplitudes  ****
    Drms = 0.0;
    count = 0;
    for(int ij=0; ij < ij_pairs; ij++) {

        if(Communicator::world->me() != ij_owner[ij])
          continue;

        i = ij_map[ij][0];
        j = ij_map[ij][1];

        for(a=0; a < pairdom_len[ij]; a++) {
          for(b=0; b < pairdom_len[ij]; b++) {
            if(fabs(T[div][ij_local[ij]][a][b]) > 1e-14 && fabs(T[dmat1][ij_local[ij]][a][b]) > 1e-14) {
              if(i != j)
                Drms += 2 * (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]) *
                       (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]);
              else
                Drms += (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]) *
                       (T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b]);
            }
            count++;
          }
        }
      //}
    }
  }
  else {
    //   ****  Compute the RMS for the New and Old Amplitudes  ****
    Drms = 0.0;
    count = 0;
    for(int ij=0; ij < ij_pairs; ij++) {
        if(Communicator::world->me() != ij_owner[ij])
          continue;

        i = ij_map[ij][0];
        j = ij_map[ij][1];
        
        for(a=0; a < pairdom_len[ij]; a++) {
          for(b=0; b < pairdom_len[ij]; b++) {
            if(fabs(T[div][ij_local[ij]][a][b]) > 1e-14 && fabs(T[dmat1][ij_local[ij]][a][b]) > 1e-14) {
              if(i != j)
                Drms += 2 * (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]) *
                       (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]);
              else
                Drms += (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]) *
                       (T_ext[nmat][ij_local[ij]][a][b] - T_ext[omat][ij_local[ij]][a][b]);
            }
            count++;
          }
        }
      //}
    }
  }

  Communicator::world->sum(Drms);
  Drms = sqrt(Drms/count);

  if(fabs(DEmp2) < econv && fabs(Drms) < rmsconv || iter >= maxiter) {
    conv = 1;
    if(iter >= maxiter)
      if(Communicator::world->me() == 0)
        fprintf(outfile, "LMP2 has not converged in the maximum number of iterations.\n maxiter = %d\n", maxiter);
  }
  else conv = 0;

}

}} // namespace psi::lmp2
