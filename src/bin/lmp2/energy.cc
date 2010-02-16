/*! \file
    \ingroup LMP2
    \brief compute the lmp2 energy
*/
//#include "mpi.h"
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <cmath>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <psifiles.h>
#include <libparallel/parallel.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;


void LMP2::energy() {
    
  int i, j, a, b, ij;
  int *ij_owner, *ij_local;
  int **ij_map;
  int **pairdomain, *pairdom_len;

  ij_map = get_ij_map();
  pairdomain = compute_pairdomain(ij_map);
  pairdom_len = compute_pairdomlen(ij_map);

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();

  if(iter > 0) 
    Emp2_old = Emp2;

  Emp2 = 0.0;
  E_SS = 0.0;
  E_OS = 0.0;

  if (iter < it_diis || diis == 0) {
    // *** Compute the new MP2 energy ***
  for(int ij=0; ij < ij_pairs; ij++) {
        if(myid != ij_owner[ij])
          continue;
        i = ij_map[ij][0];
        j = ij_map[ij][1];

        if (i!=j) {
          for(a=0; a < pairdom_len[ij]; a++) {
            for(b=0; b < pairdom_len[ij]; b++) {
              if(fabs(Ktilde[ij_local[ij]][a][b]) > 1e-14) {
                E_OS += 2 * Ktilde[ij_local[ij]][a][b] * ( T[div][ij_local[ij]][a][b] );
                E_SS += 2 * Ktilde[ij_local[ij]][a][b] * ( T[div][ij_local[ij]][a][b] - T[div][ij_local[ij]][b][a]);
                Emp2 += 2 * Ktilde[ij_local[ij]][a][b] * (2 * T[div][ij_local[ij]][a][b] - T[div][ij_local[ij]][b][a]);
              }
            }
          }
        }
        else {
          for(a=0; a < pairdom_len[ij]; a++) {
            for(b=0; b < pairdom_len[ij]; b++) {
              if(fabs(Ktilde[ij_local[ij]][a][b]) > 1e-14) {
                E_OS += 2 * Ktilde[ij_local[ij]][a][b] * ( T[div][ij_local[ij]][a][b] );
                E_SS += 2 * Ktilde[ij_local[ij]][a][b] * ( T[div][ij_local[ij]][a][b] - T[div][ij_local[ij]][b][a]);
                Emp2 += Ktilde[ij_local[ij]][a][b] * (2 * T[div][ij_local[ij]][a][b] - T[div][ij_local[ij]][b][a]);
              }
            }
          }
        }
    }
  }
  else {
    // *** Compute the new MP2 energy ***
  for(int ij=0; ij < ij_pairs; ij++) {
        if(myid != ij_owner[ij])
          continue;

        i = ij_map[ij][0];
        j = ij_map[ij][1];

        if (i!=j) {
          for(a=0; a < pairdom_len[ij]; a++) {
            for(b=0; b < pairdom_len[ij]; b++) {
              if(fabs(Ktilde[ij_local[ij]][a][b]) > 1e-14) 
                E_OS += 2 * Ktilde[ij_local[ij]][a][b] * ( T_ext[nmat][ij_local[ij]][a][b] );
                E_SS += 2 * Ktilde[ij_local[ij]][a][b] * ( T_ext[nmat][ij_local[ij]][a][b] -
                            T_ext[nmat][ij_local[ij]][b][a]);
                Emp2 += 2 * Ktilde[ij_local[ij]][a][b] * (2 * T_ext[nmat][ij_local[ij]][a][b] -
                            T_ext[nmat][ij_local[ij]][b][a]);
            }
          }
        }
        else {
          for(a=0; a < pairdom_len[ij]; a++) {
            for(b=0; b < pairdom_len[ij]; b++) {
              if(fabs(Ktilde[ij_local[ij]][a][b]) > 1e-14) 
                E_OS += 2 * Ktilde[ij_local[ij]][a][b] * ( T_ext[nmat][ij_local[ij]][a][b] );
                E_SS += 2 * Ktilde[ij_local[ij]][a][b] * ( T_ext[nmat][ij_local[ij]][a][b] -
                            T_ext[nmat][ij_local[ij]][b][a]);
                Emp2 += Ktilde[ij_local[ij]][a][b] * (2 * T_ext[nmat][ij_local[ij]][a][b] -
                            T_ext[nmat][ij_local[ij]][b][a]);
            }
          }
        }
    }
  }

  Communicator::world->sum(Emp2);

}

void LMP2::print_iteration() {

    if(iter != 0) {
      if(iter > it_diis && diis != 0) {
        fprintf(outfile, "%d\t %20.12f\t %20.12f %20.12f DIIS\n", iter, Emp2, DEmp2, Drms);
        fflush(outfile);
      }
      else  {
        fprintf(outfile, "%d\t %20.12f\t %20.12f %20.12f\n", iter, Emp2, DEmp2, Drms);
        fflush(outfile);
      }
    }
    else {
      fprintf(outfile, "\nBegin LMP2 Iterations\n");
      fprintf(outfile, "Iter\t           LMP2\t\t\t Delta(LMP2)\t\t RMS(T)\n");
      fprintf(outfile, "%d\t %20.12f\n", iter, Emp2);
    }

}

}} // namespace psi::lmp2
