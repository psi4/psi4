/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/

#include "mpi.h"
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


int LMP2::get_nso() {
  int nso_;

  if(Communicator::world->me() == 0)
    nso_ = chkpt->rd_nso();

  // This version of bcast sends a single value from master
  if(nprocs > 1)
    Communicator::world->bcast(nso_);

  return nso_;

}

int LMP2::get_natom() {
  int natom_;

  if(Communicator::world->me() == 0)
    natom_ = chkpt->rd_natom();

  // This version of bcast sends a single value from master
  if(nprocs > 1)
    Communicator::world->bcast(natom_);

  return natom_;

}

int LMP2::get_nirreps() {
  int nirreps_;

  if(Communicator::world->me() == 0)
    nirreps_ = chkpt->rd_nirreps();
  
  if(nprocs > 1)
    Communicator::world->bcast(nirreps_);

  return nirreps_;

}

int LMP2::get_nshell() {
  int nshell_;

  if(Communicator::world->me() == 0)
    nshell_ = chkpt->rd_nshell();

  if(nprocs > 1)
    Communicator::world->bcast(nshell_);

  return nshell_;

}

int LMP2::get_puream() {
  int puream_;

  if(Communicator::world->me() == 0)
    puream_ = chkpt->rd_puream();

  if(nprocs > 1)
    Communicator::world->bcast(puream_);

  return puream_;

}

int* LMP2::get_doccpi() {
  int *docc;

  if(Communicator::world->me() == 0)
    docc = chkpt->rd_clsdpi();
  else
    docc = init_int_array(nirreps);

  if(nprocs > 1)
    Communicator::world->bcast(docc, nirreps, 0);

  return docc;

}

int* LMP2::get_soccpi() {
  int *socc;

  if(Communicator::world->me() == 0)
    socc = chkpt->rd_openpi();
  else
    socc = init_int_array(nirreps);

  if(nprocs > 1)
    Communicator::world->bcast(socc, nirreps, 0);
  
  return socc;

}

int LMP2::get_frdocc() {
  int *frd;

  if(Communicator::world->me() == 0)
    frd = chkpt->rd_frzcpi();
  else
    frd = init_int_array(nirreps);

  if(nprocs > 1)
    Communicator::world->bcast(frd, nirreps, 0);
  
  return frd[0];

}


int* LMP2::get_stype() {
  int *sty;

  if(Communicator::world->me() == 0)
    sty = chkpt->rd_stype();
  else
    sty = init_int_array(nshell);

  if(nprocs > 1)
    Communicator::world->bcast(sty, nshell, 0);

  return sty;

}

int* LMP2::get_snuc() {
  int *snu;

  if(Communicator::world->me() == 0)
    snu = chkpt->rd_snuc();
  else
    snu = init_int_array(nshell);

  if(nprocs > 1)
    Communicator::world->bcast(snu, nshell, 0);

  return snu;

}

double LMP2::get_enuc() {
  double enuc_;

  if(Communicator::world->me() == 0)
    enuc_ = chkpt->rd_enuc();

  if(nprocs > 1)
    Communicator::world->bcast(enuc_);
  
  return enuc_;

}

double LMP2::get_escf() {
  double escf_;

  if(Communicator::world->me() == 0)
    escf_ = chkpt->rd_escf();

  if(nprocs > 1)
    Communicator::world->bcast(escf_);

  return escf_;

}

double* LMP2::get_evals() {

  double *evals_;

  if(Communicator::world->me() == 0)
    evals_ = chkpt->rd_evals();
  else
    evals_ = init_array(nso);

  if(nprocs > 1)
    Communicator::world->bcast(evals_, nso, 0);
  
  return evals_;

}

double** LMP2::get_MOC() {

  double **C_;

  if(Communicator::world->me() == 0)
    C_ = chkpt->rd_scf();
  else
    C_ = block_matrix(nso,nso);

  if(nprocs > 1)
    Communicator::world->bcast(C_[0], nso*nso, 0);

  return C_;

}

}}
