/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/

#include <iostream>
#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.hpp>
//#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
#include <libmints/mints.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

using namespace psi;

namespace psi{

namespace lmp2{

LMP2::LMP2() {
  fprintf(outfile, "Do not use the default constructor for LMP2");
}

LMP2::~LMP2() {
//  free_block(D);
}

LMP2::LMP2(shared_ptr<PSIO> psio_o, shared_ptr<Chkpt> chkpt_o) {

  psio = psio_o;
  chkpt = chkpt_o;

}

void LMP2::get_moinfo() {

  nso = get_nso();
  natom = get_natom();
  nirreps = get_nirreps();
  nshell = get_nshell();
  C = get_MOC();

  print_moinfo();

  return;
}

void LMP2::opdm() {

  int i;
  int *doccpi;

  doccpi = get_doccpi();

  for(i=0,nocc=0; i < nirreps; i++) nocc += doccpi[i];

  D = block_matrix(nso,nso);
  C_DGEMM('n', 't', nso, nso, nocc, 1, &(C[0][0]), nso, &(C[0][0]), nso, 0, &(D[0][0]), nso);

}

void LMP2::print_moinfo(){

  // A couple of error traps
  if(nirreps != 1) {
    if(Communicator::world->me() == 0) {
        char *symm_label = chkpt->rd_sym_label();
        throw InputException("Local MP2 is only valid in C1 symmetry", symm_label, __FILE__, __LINE__);
    }
  }


  double Enuc = get_enuc();
  Escf = get_escf();

  if(Communicator::world->me() == 0) {
    fprintf(outfile,"\n");
    fprintf(outfile,"\tChkpt Parameters:\n");
    fprintf(outfile,"\t--------------------\n");
    fprintf(outfile,"\tNumber of Irreps    \t= %d\n",nirreps);
    fprintf(outfile,"\tNumber of MOs    \t= %d\n",nso);
    fprintf(outfile,"\n");
    fprintf(outfile,"\tNuclear rep. energy   \t= %20.15f\n",Enuc);
    fprintf(outfile,"\tSCF energy            \t= %20.15f\n",Escf);
    fflush(outfile);
  }
}

void LMP2::get_fock() {

  int i, j;
  double **temp, *X;

  aoF = block_matrix(nso, nso);
  loF = block_matrix(nso, nso);
  X = init_array((nso*nso+nso)/2);
  temp = block_matrix(nso, nso);

  if(Communicator::world->me() == 0)
    X = chkpt->rd_fock();

  if(Communicator::world->nproc() > 1)
    Communicator::world->bcast(X, (nso*nso+nso)/2, 0);
  
  tri_to_sq(X, aoF, nso);
  free(X);

  // Transform the Fock matrix from the AO to LO basis
  C_DGEMM('t', 'n', nso, nso, nso, 1, &(C[0][0]), nso, &(aoF[0][0]), nso, 0, &(temp[0][0]), nso);
  C_DGEMM('n', 'n', nso, nso, nso, 1, &(temp[0][0]), nso, &(C[0][0]), nso, 0, &(loF[0][0]), nso);
  free_block(temp);

}

int* LMP2::get_ij_owner() {

  int i, j, ij, v, count;
  int *ij_owner;

  ij_owner = init_int_array(ij_pairs);

  v=0;
  for(ij=0; ij < ij_pairs; ij++) {
    ij_owner[ij] = v % Communicator::world->nproc();
    v++;
  }

  return &(ij_owner[0]);

}

/* 
** we distribute ij's among processes as follows...
**
** ij0 -> proc0
** ij1 -> proc1
** ij2 -> proc2
** ...
** ij(nprocs-1) -> proc(nprocs-1)
** ij(nprocs)   -> proc0     (start going back around)
** ij(nprocs+1) -> proc1
** ...
** ..
** This routine creates a mapping array to get us from a global ij idex
** to the corresponding local index on a particular process, so we
** know where in local memory to find it.
**
** From the above mapping, we see that 
** ij_local[0] = 0
** ij_local[1] = 0
** ..
** ij_local[nprocs-1] = 0
** ij_local[nprocs]   = 1
** ij_local[nprocs+1] = 1
** ...
*/
int* LMP2::get_ij_local() {

  int i, j, ij, v, count;
  int *ij_local;

  ij_local = init_int_array(ij_pairs);

  v=0;
  count=0;
  for(ij=0; ij < ij_pairs; ij++) {
      ij_local[ij] = count;
      if(v%Communicator::world->nproc() == Communicator::world->nproc()-1) count++;
      v++;
  }

  return &(ij_local[0]);

}

int* LMP2::get_mn_owner(int n) {

  int count, v, num_unique_shells;
  int *mn_owner_;

  num_unique_shells = n;

  mn_owner_ = init_int_array(n);

  v = 0;
  for(count=0; count < num_unique_shells; count++) {
    mn_owner_[count] = v%Communicator::world->nproc();
    v++;
  }

  return &(mn_owner_[0]);

}

int LMP2::get_mn_pairs(int n) {

  int count, v, num_unique_shells;
  int mn_pairs_;

  num_unique_shells = n;

  v = 0;
  for(count=0; count < num_unique_shells; count++) {
    if(Communicator::world->me()== v% Communicator::world->nproc()) mn_pairs_++;
    v++;
  }

  return mn_pairs_;
}

int LMP2::get_num_unique_shells() {

    int count = 0;
    for (int M = 0; M < nshell; M++) {
        for (int N = 0; N <= M; N++, count++) {
        }
    }

    return count;
}

int** LMP2::get_MN_shell(shared_ptr<BasisSet> basisset) {

    int ** MN_shell = init_int_matrix(4, num_unique_shells);
    int count = 0;
    for (int M = 0; M < nshell; M++) {
        int numm = basisset->shell(M)->nfunction();
        for (int N = 0; N <= M; N++, count++) {
            int numn = basisset->shell(N)->nfunction();
            MN_shell[0][count] = M;
            MN_shell[1][count] = numm;
            MN_shell[2][count] = N;
            MN_shell[3][count] = numn;
        }
    }

    return MN_shell;
}

int** LMP2::get_ij_map() {

    int i, j, ij;
    int **ij_map_ = init_int_matrix(ij_pairs,2);
    int counter;


    counter = 0;
    for (i = 0, ij=0 ; i < nocc; i++) {
        for (j = 0; j <= i; j++, ij++) {
            if (pairdom_exist[ij]) {
                ij_map_[counter][0] = i;
                ij_map_[counter][1] = j;
                counter++;
            }
        }
    }

    //for (int ij2 = 0; ij2<ij_pairs; ij2++)
	//fprintf(outfile," %d (%d %d)\n", ij2,ij_map_[ij2][0], ij_map_[ij2][1]);

    return ij_map_;

}

// don't try to use this for ij's which have been neglected
int* LMP2::original_ij_map() {

    int pairs = (nocc * (nocc + 1)) / 2;
    int i, j, ij;
    int *map_ = init_int_array(pairs);
    int counter;


    counter = 0;
    for (i = 0, ij=0 ; i < nocc; i++) {
        for (j = 0; j <= i; j++, ij++) {
            if (pairdom_exist[ij]) {
            	map_[ij] = counter;
                counter++;
            }
	    else 
		map_[ij] = -1;
        }
    }

    return map_;

}

int **LMP2::compute_pairdomain(int **ij_map_) {

    int m, i, j, ij, a;
    //int **ij_map_ = get_ij_map(&(pairdom_exist_[0]));
    int **pairdomain_ = init_int_matrix(ij_pairs, natom);

    for(m=0; m < ij_pairs; m++) {
        i = ij_map_[m][0];
        j = ij_map_[m][1];
        ij = (i * (i + 1)) / 2 + j;
        if (pairdom_exist[ij]) {
            for (a=0; a < natom; a++) {
                if (domain[i][a] || domain[j][a]) {
                    pairdomain_[m][a] = 1;
                }
            }
        }
    }

    return pairdomain_;
}

int *LMP2::compute_pairdomlen(int **ij_map_) {

    int m, i, j, ij, a;
    //int **ij_map_ = get_ij_map(&(pairdom_exist_[0]));
    int *pairdom_len_ = init_int_array(ij_pairs);

    for(m=0; m < ij_pairs; m++) {
        i = ij_map_[m][0];
        j = ij_map_[m][1];
        ij = (i * (i + 1)) / 2 + j;
        if (pairdom_exist[ij]) {
            for (a=0; a < natom; a++) {
                if (domain[i][a] || domain[j][a]) {
                    pairdom_len_[m] += aostop[a] - aostart[a] + 1;
                }
            }
        }
    }

    return pairdom_len_;
}

}} /* End namespace */


