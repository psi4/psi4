#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdio>

#include <liboptions/liboptions.h>
#include <libutil/libutil.h>

#include "scf.h"

extern FILE* outfile;

using namespace boost;

namespace psi{
namespace mcscf{

extern MemoryManager* memory_manager;

SCF::SCF(Options& options_, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
: Wavefunction(options_, psio, chkpt)
{
//  startup();
}

SCF::~SCF()
{
    // This is now called explicity after a run, because the pointer to the wavefunction object
    // is held in Process::environment, so this doesn't go out of scope until the end of the run
  //cleanup();
}

void SCF::startup()
{
  ioff    = moinfo_scf->get_ioff();
  nirreps = moinfo_scf->get_nirreps();
  nso     = moinfo_scf->get_nso();
  sopi    = moinfo_scf->get_sopi();

  docc    = moinfo_scf->get_docc();
  actv    = moinfo_scf->get_actv();

  // Compute the block_offset for the so basis
  allocate1(double,block_offset,nirreps);
  block_offset[0] = 0;
  for(int h=1; h < nirreps; ++h){
    block_offset[h] = block_offset[h-1] + sopi[h-1];
  }

  // Allocate the pairs
  allocate1(int,pairpi,nirreps);
  allocate1(int,pair_offset,nirreps);
  allocate2(int,pair,nso,nso);
  allocate2(int,pair_sym,nso,nso);
  pairs = 0;
  PK = 0;
  K = 0;

  if(options_.get_str("REFERENCE")=="RHF"){
    reference = rhf;
    same_orbs_ = true;
    same_dens_ = true;
  }
  else if(options_.get_str("REFERENCE")=="ROHF"){
    reference = rohf;
    same_orbs_ = true;
    same_dens_ = false;
  }
  else if(options_.get_str("REFERENCE")=="UHF"){
    reference = uhf;
    same_orbs_ = false;
    same_dens_ = false;
  }
  else if(options_.get_str("REFERENCE")=="TWOCON"){
    reference = tcscf;
    same_orbs_ = true;
    same_dens_ = false;
    if(moinfo_scf->get_guess_occupation()){
      printf("\n  ERROR:  MCSCF cannot guess the active orbital occupation\n");
      fprintf(outfile,"\n\n  MCSCF cannot guess the active orbital occupation\n");
      fflush(outfile);
      exit(1);
    }
  }

  root = options_.get_int("FOLLOW_ROOT") - 1;

  // OUT OF CORE ALGORITHM
  out_of_core  = false;

  // DIIS
  use_diis = options_.get_bool("USE_DIIS");
  ndiis    = options_.get_int("NDIIS");
  current_diis = 0;

  turn_on_actv = options_.get_int("TURN_ON_ACTV");

  epsilon     .allocate("epsilon",nirreps,sopi);

  e           .allocate("e",nirreps,sopi,sopi);
  C           .allocate("C",nirreps,sopi,sopi);
  C_t         .allocate("C_t",nirreps,sopi,sopi);
  C_T         .allocate("C_T",nirreps,sopi,sopi);
  Dc          .allocate("Dc",nirreps,sopi,sopi);
  Dc_old          .allocate("Dc_old",nirreps,sopi,sopi);
  Feff_t      .allocate("Feff_t",nirreps,sopi,sopi);
  Feff_t_old  .allocate("Feff_t",nirreps,sopi,sopi);
  Feff_oAO    .allocate("Feff_oAO",nirreps,sopi,sopi);
  Fc          .allocate("Fc",nirreps,sopi,sopi);
  Fc_t        .allocate("Fc_t",nirreps,sopi,sopi);
  G           .allocate("G",nirreps,sopi,sopi);
  H           .allocate("H",nirreps,sopi,sopi);
  O           .allocate("O",nirreps,sopi,sopi);
  S           .allocate("S",nirreps,sopi,sopi);
  S_sqrt_inv  .allocate("S^-1/2",nirreps,sopi,sopi);
  S_sqrt      .allocate("S^1/2",nirreps,sopi,sopi);
  T           .allocate("T",nirreps,sopi,sopi);

  for(int i = 0; i < ndiis; ++i){
    diis_F[i] .allocate("diis_F[" + to_string(i) + "]",nirreps,sopi,sopi);
    diis_e[i] .allocate("diis_e[" + to_string(i) + "]",nirreps,sopi,sopi);
  }

  if(reference == rohf){
    Do        .allocate("Do",nirreps,sopi,sopi);
    Fo        .allocate("Fo",nirreps,sopi,sopi);
    Fo_t      .allocate("Fo_t",nirreps,sopi,sopi);
  }
  if(reference == tcscf){
    int count = 0;
    for(int h = 0; h < nirreps; ++h){
      for(int n = 0; n < actv[h]; ++n){
        tcscf_sym[count] = h;
        tcscf_mos[count] = docc[h] + n;
        count++;
      }
    }

    nci = count;
    fprintf(outfile,"\n  TWOCON MOs = [");
    for(int I = 0; I < nci; ++I)
      fprintf(outfile,"%d (%s)%s",tcscf_mos[I] + block_offset[tcscf_sym[I]],
                                 moinfo_scf->get_irr_labs(tcscf_sym[I]),
                                 I != nci - 1 ? "," : "");
    fprintf(outfile,"]");

    Favg      .allocate("Favg",nirreps,sopi,sopi);
    Favg_t    .allocate("Favg_t",nirreps,sopi,sopi);

    allocate1(double,ci,nci);
    allocate1(double,ci_grad,nci);
    allocate2(double,H_tcscf,nci,nci);
    for(int I = 0; I < nci; ++I){
      Dtc[I]    .allocate("Dtc[" + to_string(I) + "]",nirreps,sopi,sopi);
      Dtc_old[I]    .allocate("Dtc_old[" + to_string(I) + "]",nirreps,sopi,sopi);
      Dsum[I]   .allocate("Dsum[" + to_string(I) + "]",nirreps,sopi,sopi);
      Ftc[I]    .allocate("Ftc[" + to_string(I) + "]",nirreps,sopi,sopi);
      Ftc_t[I]  .allocate("Ftc_t[" + to_string(I) + "]",nirreps,sopi,sopi);
      ci[I] = 0.0;// (I == 0 ? 0.7071067811865475244 : -0.7071067811865475244) ;
    }
    if(options_.get_bool("FORCE_TWOCON")){
      ci[0] =   0.7071067811865475244;
      ci[1] = - 0.7071067811865475244;
    }
  }
}

void SCF::cleanup()
{
  release1(block_offset);
  release1(pairpi);
  release1(pair_offset);
  release1(pairs);
  release2(pair);
  release2(pair_sym);

  if(epsilon.is_allocated())    epsilon.subtract_reference();
  if(C.is_allocated())          C.subtract_reference();
  if(C_t.is_allocated())        C_t.subtract_reference();
  if(C_T.is_allocated())        C_T.subtract_reference();
  if(Dc.is_allocated())         Dc.subtract_reference();
  if(Dc_old.is_allocated())     Dc_old.subtract_reference();
  if(Do.is_allocated())         Do.subtract_reference();
  if(Fc.is_allocated())         Fc.subtract_reference();
  if(Fc_t.is_allocated())       Fc_t.subtract_reference();
  if(Fo.is_allocated())         Fo.subtract_reference();
  if(Fo_t.is_allocated())       Fo_t.subtract_reference();
  if(Favg.is_allocated())       Favg.subtract_reference();
  if(Favg_t.is_allocated())     Favg_t.subtract_reference();
  if(Feff_t.is_allocated())     Feff_t.subtract_reference();
  if(Feff_t_old.is_allocated()) Feff_t_old.subtract_reference();
  if(Feff_oAO.is_allocated())   Feff_oAO.subtract_reference();
  if(G.is_allocated())          G.subtract_reference();
  if(T.is_allocated())          T.subtract_reference();
  if(H.is_allocated())          H.subtract_reference();
  if(O.is_allocated())          O.subtract_reference();
  if(S.is_allocated())          S.subtract_reference();
  if(S_sqrt_inv.is_allocated()) S_sqrt_inv.subtract_reference();
  if(S_sqrt.is_allocated())     S_sqrt.subtract_reference();
  if(e.is_allocated())          e.subtract_reference();
  for(int i = 0; i < maxci; ++i){
      if(Ftc[i].is_allocated())   Ftc[i].subtract_reference();
      if(Ftc_t[i].is_allocated()) Ftc_t[i].subtract_reference();
      if(Dtc[i].is_allocated())   Dtc[i].subtract_reference();
      if(Dtc_old[i].is_allocated())   Dtc_old[i].subtract_reference();
      if(Dsum[i].is_allocated())  Dsum[i].subtract_reference();
  }
  for(int i = 0; i < maxdiis; ++i){
      if(diis_F[i].is_allocated()) diis_F[i].subtract_reference();
      if(diis_e[i].is_allocated()) diis_e[i].subtract_reference();
  }

  if(reference == tcscf){
    release1(ci);
    release1(ci_grad);
    release2(H_tcscf);
  }

  release1(PK);
  release1(K);
}

void SCF::transform(SBlockMatrix& Initial, SBlockMatrix& Final, SBlockMatrix& Transformation){
  T.multiply(false,false,Initial,Transformation);
  Final.multiply(true,false,Transformation,T);
}

}} /* End Namespaces */
