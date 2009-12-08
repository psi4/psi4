/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/

#include "mpi.h"
#include <iostream>
#include <fstream>              // file I/O support
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

int LMP2::get_nso() {
  int nso_;

  if(myid == 0) {
    nso_ = chkpt->rd_nso();
    MPI_Bcast(&nso_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Bcast(&nso_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return nso_;

}

int LMP2::get_natom() {
  int natom_;

  if(myid == 0) {
    natom_ = chkpt->rd_natom();
    MPI_Bcast(&natom_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Bcast(&natom_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return natom_;

}

int LMP2::get_nirreps() {
  int nirreps_;

  if(myid == 0) {
    nirreps_ = chkpt->rd_nirreps();
    MPI_Bcast(&nirreps_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Bcast(&nirreps_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return nirreps_;

}

int LMP2::get_nshell() {
  int nshell_;

  if(myid == 0) {
    nshell_ = chkpt->rd_nshell();
    MPI_Bcast(&nshell_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else 
    MPI_Bcast(&nshell_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return nshell_;

}

int LMP2::get_puream() {
  int puream_;

  if(myid == 0) {
    puream_ = chkpt->rd_puream();
    MPI_Bcast(&puream_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else 
    MPI_Bcast(&puream_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return puream_;

}

int* LMP2::get_doccpi() {
  int *docc;

  if(myid == 0) {
    docc = chkpt->rd_clsdpi();
    MPI_Bcast(&(docc[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    docc = init_int_array(nirreps);
    MPI_Bcast(&(docc[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return docc;

}

int* LMP2::get_soccpi() {
  int *socc;

  if(myid == 0) {
    socc = chkpt->rd_openpi();
    MPI_Bcast(&(socc[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    socc = init_int_array(nirreps);
    MPI_Bcast(&(socc[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return socc;

}

int LMP2::get_frdocc() {
  int *frd;

  if(myid == 0) {
    frd = chkpt->rd_openpi();
    MPI_Bcast(&(frd[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    frd = init_int_array(nirreps);
    MPI_Bcast(&(frd[0]), nirreps, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return frd[0];

}


int* LMP2::get_stype() {
  int *sty;

  if(myid == 0) {
    sty = chkpt->rd_stype();
    MPI_Bcast(&(sty[0]), nshell, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else  {
    sty = init_int_array(nshell);
    MPI_Bcast(&(sty[0]), nshell, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return sty;

}

int* LMP2::get_snuc() {
  int *snu;

  if(myid == 0) {
    snu = chkpt->rd_snuc();
    MPI_Bcast(&(snu[0]), nshell, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    snu = init_int_array(nshell);
    MPI_Bcast(&(snu[0]), nshell, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return snu;

}

double LMP2::get_enuc() {
  double enuc_;

  if(myid == 0) {
    enuc_ = chkpt->rd_enuc();
    MPI_Bcast(&enuc_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else 
    MPI_Bcast(&enuc_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return enuc_;

}

double LMP2::get_escf() {
  double escf_;

  if(myid == 0) {
    escf_ = chkpt->rd_escf();
    MPI_Bcast(&escf_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else 
    MPI_Bcast(&escf_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return escf_;

}

double* LMP2::get_evals() {

  double *evals_;

  if(myid == 0) {
    evals_ = chkpt->rd_evals();
    MPI_Bcast(&(evals_[0]), nso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    evals_ = init_array(nso);
    MPI_Bcast(&(evals_[0]), nso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  return evals_;

}

double** LMP2::get_MOC() {

  double **C_;

  if(myid == 0) {
    C_ = chkpt->rd_scf();
    MPI_Bcast(&(C_[0][0]), nso*nso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    C_ = block_matrix(nso,nso);
    MPI_Bcast(&(C_[0][0]), nso*nso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  return C_;

}

}}
