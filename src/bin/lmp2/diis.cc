/*! \file
    \ingroup LMP2
    \brief compute the lmp2 energy
*/
#include "mpi.h"
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;

    
void LMP2::diis_ext() {

  int i, j, k, p, q, a, b, v, ij;
  double **Bmat, **co;
  int *work;
  int c_ext;
  double ***temp;
  int *ij_local, *ij_owner;
  int **ij_map;
  int **pairdomain, *pairdom_len;

  ij_map = get_ij_map();
  pairdomain = compute_pairdomain(ij_map);
  pairdom_len = compute_pairdomlen(ij_map);

  temp = (double ***) malloc(ij_pairs * sizeof(double **));
  for(ij=0; ij < ij_pairs; ij++) {
    temp[ij] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
  }

  ij_local = get_ij_local();
  ij_owner = get_ij_owner();

    for (ij = 0; ij < ij_pairs; ij++) {
      if(myid != ij_owner[ij])
        continue;
      for(a=0; a < pairdom_len[ij]; a++) {
        for(b=0; b < pairdom_len[ij]; b++) {
          temp[ij][a][b] = T[div][ij_local[ij]][a][b] - T[dmat1][ij_local[ij]][a][b];
        }
      }
    //}
  }

  for(ij=0; ij < ij_pairs; ij++) {
    Communicator::world->sum(temp[ij][0], pairdom_len[ij]*pairdom_len[ij], error[dmat2][ij][0]);
  }

  for(ij=0; ij < ij_pairs; ij++) {
    free_block(temp[ij]);
  }
  free(temp);
/*  if(myid == 0) {
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++) {
      fprintf(outfile, "error[%d][%d]:\n", dmat2, ij);
      print_mat(error[dmat2][ij], pairdom_len[ij], pairdom_len[ij], outfile);
  }}
  }
*/


  if(iter >= it_diis) {

    // *** Compute the B matrix ***
    Bmat = block_matrix(matsize+1,matsize+1);
    for(p=0; p <= matsize; p++) {
      for(q=0; q <= matsize; q++) {
        if(p==matsize || q==matsize) 
          Bmat[p][q] = -1.0;
        else {
          Bmat[p][q] = 0.0;
            for (ij = 0; ij < ij_pairs; ij++) {
              for(a=0; a < pairdom_len[ij]; a++) {
                for(b=0; b < pairdom_len[ij]; b++) {
                  Bmat[p][q] += error[p][ij][a][b] * error[q][ij][a][b];
                }
              }
            //}
          }
        }
      }
    }
    Bmat[matsize][matsize] = 0.0;
/*    if(myid == 0) {
      fprintf(outfile, "Bmat:\n");
      print_mat(Bmat, matsize+1, matsize+1, outfile);
    }
*/

    co = block_matrix(matsize+1,1);
    co[matsize][0] = -1.0;

    work = init_int_array(matsize+1);

                                                  // use DGESV to compute the solution to Ax=B
    C_DGESV(matsize+1, 1, Bmat[0], matsize+1, work, co[0], matsize+1);

                                                     // compute the new MP2 coefficients
    v=0;
    for (ij = 0; ij < ij_pairs; ij++, v++) {
        if(v%nprocs != myid)
          continue;
        for(a=0; a < pairdom_len[ij]; a++) {
          for(b=0; b < pairdom_len[ij]; b++) {
            T_ext[nmat][ij_local[ij]][a][b] = 0.0;
            for(k=0; k < matsize; k++) {
              if(k == ndiis - 1) c_ext = k-k + 1;
              else if(k == ndiis - 2) c_ext = k - k;
              else c_ext = k + 2;
                T_ext[nmat][ij_local[ij]][a][b] += co[k][0] * T[c_ext][ij_local[ij]][a][b];
            }
          }
        }
      //}
    }

/*    v=0;
    for(i=0, ij=0; i < nocc; i++) {
      for(j=0; j <= i; j++, ij++, v++) {
        if(v%nprocs != myid)
          continue;
        fprintf(outfile, "Tdiis[%d] Matrix :\n", ij);
        print_mat(T_ext[nmat][ij], pairdom_len[ij], pairdom_len[ij], outfile);
    }} 
*/

    free_block(Bmat);
    free_block(co);
    free(work);
  }


}

}} // namespace psi::lmp2
