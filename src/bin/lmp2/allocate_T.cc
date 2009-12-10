/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdlib>              // support for exit()
//#include <cstring>
//#include <memory.h>
//#include <assert.h>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;

void LMP2::allocate_T() {



  int i, j, ij;
  int *ij_local, *ij_owner;

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();

  if(diis == 1) {
    // Allocate memory for the amplitudes
    T = (double ****) malloc(ndiis * sizeof(double ***));
    for(i=0; i < ndiis; i++) {
      if(ij_pairs%nprocs == 0) {
        T[i] = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
        for(ij=0; ij < ij_pairs; ij++) {
          if(myid == ij_owner[ij])
            T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        }
      }
      else {
        if(myid < ij_pairs%nprocs) {
          T[i] = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
        else {
          T[i] = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
      }
    }

    // Allocate memory for the error matrices
    error = (double ****) malloc(ndiis * sizeof(double ***));
    for(i=0; i < ndiis; i++) {
      error[i] = (double ***) malloc(ij_pairs * sizeof(double **));
      for(j=0; j < ij_pairs; j++) {
        error[i][j] = block_matrix(pairdom_len[j], pairdom_len[j]);
      }
    }

    // Allocate memory for the DIIS extrapolated amplitudes
    T_ext = (double ****) malloc(2 * sizeof(double ***));
    for(i=0; i < 2; i++) {

      if(ij_pairs%nprocs == 0) {
        T_ext[i] = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
        for(ij=0; ij < ij_pairs; ij++) {
          if(myid == ij_owner[ij])
            T_ext[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        }
      }
      else {
        if(myid < ij_pairs%nprocs) {
          T_ext[i] = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T_ext[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
        else {
          T_ext[i] = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T_ext[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
      }

    }
  }
  else {
    // Allocate memory for the amplitudes
    T = (double ****) malloc(2 * sizeof(double ***));
    for(i=0; i < 2; i++) {

      if(ij_pairs%nprocs == 0) {
        T[i] = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
        for(ij=0; ij < ij_pairs; ij++) {
          if(myid == ij_owner[ij])
            T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        }
      }
      else {
        if(myid < ij_pairs%nprocs) {
          T[i] = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
        else {
          T[i] = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
          for(ij=0; ij < ij_pairs; ij++) {
            if(myid == ij_owner[ij])
              T[i][ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
          }
        }
      }
    }
  }


}

}} // namespace psi::lmp2
