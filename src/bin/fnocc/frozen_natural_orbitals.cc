/**
  * Frozen natural orbitals
  * Eugene DePrince
  * June 2012
  *
  */

#include"psi4-dec.h"
#include<psifiles.h>
#include<libmints/mints.h>
#include<libmints/mintshelper.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libtrans/mospace.h>
#include<libtrans/integraltransform.h>
#include<libiwl/iwl.h>
#include"ccsd.h"
#include"blas.h"
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>
#include<lib3index/dftensor.h>
#include<lib3index/cholesky.h>


using namespace psi;
using namespace boost;

namespace psi{namespace fnocc{
void SortOVOV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt);
}}

namespace psi{namespace fnocc{

FrozenNO::FrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options):
  Wavefunction(options, _default_psio_lib_)
{
    // copy wave function.
    copy(wfn);
    common_init();
}
FrozenNO::~FrozenNO()
{
}

void FrozenNO::common_init() {
    if (nirrep_>1){
       throw PsiException("FrozenNO requires symmetry c1 (for now!)",__FILE__,__LINE__);
    }
    nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
    for (int h=0; h<nirrep_; h++){
        nfzc   += frzcpi_[h];
        nfzv   += frzvpi_[h];
        nso    += nsopi_[h];
        nmo    += nmopi_[h];
        ndocc  += doccpi_[h];
    }
    ndoccact = ndocc - nfzc;
    nvirt    = nmo - ndocc;

    // quit if number of virtuals is less than number of doubly occupied
    if (nvirt<ndoccact){
       throw PsiException("ndocc must be less than nvirt",__FILE__,__LINE__);
    }

}
// use this function to return the mp2 energy in the full basis.
double FrozenNO::compute_energy(){
  return emp2;
}

/*
 * build natural orbitals and transform TEIs
 */
void FrozenNO::ComputeNaturalOrbitals(){

  long int o = ndoccact;
  long int v = nvirt;

  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *               Frozen Natural Orbitals               *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  // transform (ia|jb) integrals
  tstart();
  TransformOVOV();
  tstop();

  // orbital energies
  tstart();
  fprintf(outfile,"        ==> Build MP2 amplitudes, OPDM, and NOs <==\n");
  fprintf(outfile,"\n");
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();
  boost::shared_ptr<Vector> eps_test = ref->epsilon_a();
  double*tempeps = eps_test->pointer();
  double *F  = tempeps + nfzc;

  long int memory = Process::environment.get_memory();
  if ( memory < 16L*o*o*v*v ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  // allocate memory for a couple of buffers 
  double*amps1 = (double*)malloc(o*o*v*v*sizeof(double));
  double*amps2 = (double*)malloc(o*o*v*v*sizeof(double));

  // build mp2 amplitudes for mp2 density
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&amps2[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  long int ijab = 0;
  emp2 = 0.0;
  double emp2_os = 0.0;
  double emp2_ss = 0.0;
  for (int a=o; a<o+v; a++){
      double da = F[a];
      for (int b=o; b<o+v; b++){
          double dab = da + F[b];
          for (int i=0; i<o; i++){
              double dabi = dab - F[i];
              for (int j=0; j<o; j++){
                  int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  double dijab = dabi-F[j];
                  amps1[ijab++] = - amps2[iajb]/dijab;
                  emp2_os -= amps2[iajb] * (amps2[iajb])/dijab;
                  emp2_ss -= amps2[iajb] * (amps2[iajb] - amps2[j*o*v*v+(a-o)*o*v+i*v+(b-o)])/dijab;
              }
          }
      }
  }
  emp2 = emp2_os + emp2_ss;

  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + Process::environment.globals["SCF TOTAL ENERGY"];

  ijab = 0;
  for (int a=o; a<o+v; a++){
      for (int b=o; b<o+v; b++){
          for (int i=0; i<o; i++){
              for (int j=0; j<o; j++){
                  int ijba = (b-o)*o*o*v+(a-o)*o*o+i*o+j;
                  amps2[ijab] = 2.0*amps1[ijab] - amps1[ijba];
                  ijab++;
              }
          }
      }
  }

  // build ab block of the density:
  double*Dab = (double*)malloc(v*v*sizeof(double));
  F_DGEMM('t','n',v,v,v*o*o,2.0,amps1,v*o*o,amps2,v*o*o,0.0,Dab,v);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);
  // reorder transformation matrix:
  double*temp    = (double*)malloc(nso*v*sizeof(double));
  for (int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("OCC_TOLERANCE");

  nvirt_no = 0;
  for (int i=0; i<v; i++) if (eigvalDab[i]>cutoff) nvirt_no++;

  fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
  fprintf(outfile,"\n");
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5li\n",v);
  fprintf(outfile,"        Number of virtual orbitals in truncated space: %5li\n",nvirt_no);
  fprintf(outfile,"\n");

  // transform Fock matrix to MP2 NO basis
  double*newFock = (double*)malloc(v*v*sizeof(double));
  memset((void*)newFock,'\0',v*v*sizeof(double));
  F_DCOPY(v,F+ndoccact,1,newFock,v+1);
  F_DGEMM('n','n',v,nvirt_no,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',nvirt_no,nvirt_no,v,1.0,temp,v,Dab,v,0.0,newFock,nvirt_no);

  // diagonalize new Fock matrix for semi-canonical orbitals
  double*neweps = (double*)malloc(nvirt_no*sizeof(double));
  Diagonalize(nvirt_no,newFock,neweps);
  // construct full mo -> no transformation matrix
  F_DGEMM('n','n',v,nvirt_no,nvirt_no,1.0,temp,v,newFock,nvirt_no,0.0,Dab,v);

  // put orbital energies back in F
  F_DCOPY(nvirt_no,neweps,1,F+ndoccact,1);

  // free memory before using libtrans
  free(temp);
  free(neweps);
  free(amps1);
  free(amps2);
  free(eigvalDab);
  free(newFock);
  tstop();

  // so->no integral transformation
  tstart();
  TransformIntegrals(Dab);
  tstop();

  // change number of frozen virtual orbitals
  frzvpi_[0] = nfzv+(nvirt-nvirt_no);

  free(Dab);
}
void FrozenNO::TransformIntegrals(double*Dab){

  long int v = nvirt;

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();

  // build our own C matrices so librans can transform so->no where sizeof(no) != sizeof(mo)
  fprintf(outfile,"        ==> Transform all two-electron integrals <==\n");
  fprintf(outfile,"\n");
  SharedMatrix Cc (new Matrix("FROZEN_OCC",nso,nfzc));
  SharedMatrix Ci (new Matrix("ACTIVE_OCC",nso,ndoccact));
  SharedMatrix Ca (new Matrix("ACTIVE_VIR",nso,nvirt_no));
  SharedMatrix Cv (new Matrix("FROZEN_VIR",nso,0));

  boost::shared_ptr<Matrix> Caomo = ref->Ca();
  double**Capointer = Caomo->pointer();

  double**pointer;
  pointer = Cc->pointer();
  for (int i=0; i<nso; i++){
      for (int j=0; j<nfzc; j++){
          pointer[i][j] = Capointer[i][j];
      }
  }
  pointer = Ci->pointer();
  for (int i=0; i<nso; i++){
      for (int j=0; j<ndoccact; j++){
          pointer[i][j] = Capointer[i][j+nfzc];
      }
  }

// so->no
  double*temp = (double*)malloc(nso*nvirt_no*sizeof(double));
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<nvirt_no; j++){
          double dum = 0.0;
          for (long int k=0; k<v; k++){
              dum += Capointer[i][ndocc+k] * Dab[j*v+k];
          }
          temp[i*nvirt_no+j] = dum;
      }
  }
  pointer = Ca->pointer();
  for (int i=0; i<nso; i++){
      for (int j=0; j<nvirt_no; j++){
          pointer[i][j] = temp[i*nvirt_no+j];
      }
  }
  free(temp);
  //pointer = Cv->pointer();
  //for (int i=0; i<nso; i++){
  //    for (int j=0; j<nfzv+(nvirt-nvirt_no); j++){
  //        pointer[i][j] = Ca[i][j+ndocc+nvirt];
  //    }
  //}

  // so->no integral transformation
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  IntegralTransform ints2(Cc, Ci, Ca, Cv, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::PitzerOrder, IntegralTransform::OccAndVir, true);
  ints2.set_keep_dpd_so_ints(1);
  ints2.set_keep_iwl_so_ints(1);
  ints2.initialize();
  ints2.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  fprintf(outfile,"\n");
}

void FrozenNO::TransformOVOV(){
  fprintf(outfile,"        ==> Transform (OV|OV) integrals <==\n");
  fprintf(outfile,"\n");
  boost::shared_ptr<psi::Wavefunction> wfn = Process::environment.wavefunction();
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::occ);
  spaces.push_back(MOSpace::vir);
  boost::shared_ptr<IntegralTransform>
      ints(new IntegralTransform(wfn,
                                 spaces,
                                 IntegralTransform::Restricted,
                                 IntegralTransform::IWLOnly,
                                 IntegralTransform::QTOrder,
                                 IntegralTransform::None,
                                 false));
  ints->set_keep_dpd_so_ints(1);
  ints->set_keep_iwl_so_ints(1);
  ints->initialize();
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
  fprintf(outfile,"\n");

  // sort integrals
  struct iwlbuf Buf;
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
  SortOVOV(&Buf,nfzc,nfzv,nmo,ndoccact,nvirt);
  iwl_buf_close(&Buf,1);
}

// DF FNO class members

DFFrozenNO::DFFrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options):
  FrozenNO(wfn,options)
{
}
DFFrozenNO::~DFFrozenNO()
{
}

// use this function to return the mp2 energy in the full basis.  
// this isn't used anymore ... i'm reluctant to remove it, though
double DFFrozenNO::compute_energy(){
  return emp2;
}

/*
 * build natural orbitals and transform TEIs
 */
void DFFrozenNO::ComputeNaturalOrbitals(){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *               Frozen Natural Orbitals               *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  fprintf(outfile,"\n");
  fprintf(outfile,"  Build 3-index integrals:\n");
  fprintf(outfile,"\n");

  long int o = ndoccact;
  long int v = nvirt;
  long int nQ;

  if ( ( options_.get_str("DF_BASIS_CC") != "CHOLESKY" ) ){
      boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
      boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, molecule(), "DF_BASIS_CC");

      boost::shared_ptr<DFTensor> DF (new DFTensor(basisset(),auxiliary,Ca(),ndocc,nvirt+nfzv,ndoccact,nvirt,options_));
      nQ = auxiliary->nbf();
      boost::shared_ptr<Matrix> tmp = DF->Qso();
      double ** Qso = tmp->pointer();

      // write Qso to disk
      boost::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_QSO,"qso",(char*)&Qso[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
  }else{
      // Cholesky 3-index integrals
      boost::shared_ptr<BasisSet> primary = basisset();
      boost::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
      double tol = options_.get_double("CHOLESKY_TOLERANCE");
      boost::shared_ptr<CholeskyERI> Ch (new CholeskyERI(boost::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0,tol,Process::environment.get_memory()));
      Ch->choleskify();
      nQ  = Ch->Q();
      boost::shared_ptr<Matrix> L = Ch->L();
      double ** Lp = L->pointer();

      // write Qso to disk
      boost::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_QSO,"qso",(char*)&Lp[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
  }

  long int memory = Process::environment.get_memory();
  if ( memory < 8L*(3L*nso*nso+nso*nso*nQ+o*v*nQ) ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  double * tmp2 = (double*)malloc(nso*nso*nQ*sizeof(double));
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_QSO,"qso",(char*)&tmp2[0],nQ*nso*nso*sizeof(double));
  psio->close(PSIF_DCC_QSO,1);

  // transform 3-index integrals from AO to MO and build Fock matrix in MO basis
  double * dfFock = (double*)malloc(nso*nso*sizeof(double));
  BuildFock(nQ,tmp2,dfFock);
  double * Qov = (double*)malloc(o*v*nQ*sizeof(double));
  for (long int q = 0; q < nQ; q++) {
      for (long int i = 0; i < o; i++) {
          for (long int a = 0; a < v; a++) {
              Qov[q*o*v+i*v+a] = tmp2[q*nmo*nmo+(i+nfzc)*nmo+(a+ndocc)];
          }
      }
  }

  // the resulting Fock matrix is NOT diagonal in the o-o and v-v spaces.  make it so!

  // virtual space:
  nvirt_no = nvirt;
  double*neweps = (double*)malloc(nvirt_no*sizeof(double));
  double*newFock = (double*)malloc(v*v*sizeof(double));
  memset((void*)newFock,'\0',v*v*sizeof(double));
  for (long int a = 0; a < v; a++) {
      for (long int b = 0; b < v; b++) {
          newFock[a*v+b] = dfFock[(a+ndocc)*nmo+(b+ndocc)];
      }
  }
  Diagonalize(nvirt_no,newFock,neweps);

  // put orbital energies back in F
  boost::shared_ptr<Vector> eps_test = epsilon_a();
  double * tempeps = eps_test->pointer();
  double * F       = tempeps + nfzc;
  double * temp    = (double*)malloc(nso*v*sizeof(double));
  double * Dab     = (double*)malloc(v*v*sizeof(double));
  F_DCOPY(nvirt_no,neweps,1,F+ndoccact,1);

  // modify c matrix
  ModifyCa(newFock);

  // transform Qov:
  for (long int q = 0; q < nQ; q++) {
     for (long int i = 0; i < o; i++) {
         for (long int a = 0; a < v; a++) {
             double dum = 0.0;
             for (long int b = 0; b < v; b++) {
                 dum += Qov[q*o*v+i*v+b] * newFock[a*v+b];
             }
             tmp2[q*o*v+i*v+a] = dum;
         }
     }
  }

  // occupied space:
  memset((void*)newFock,'\0',o*o*sizeof(double));
  for (long int i = 0; i < o; i++) {
      for (long int j = 0; j < o; j++) {
          newFock[i*o+j] = dfFock[(i+nfzc)*nmo+(j+nfzc)];
      }
  }
  Diagonalize(o,newFock,neweps);

  // put orbital energies back in F
  F_DCOPY(o,neweps,1,F,1);

  // modify c matrix
  ModifyCa_occ(newFock);

  // transform Qov:
  for (long int q = 0; q < nQ; q++) {
     for (long int i = 0; i < o; i++) {
         for (long int a = 0; a < v; a++) {
             double dum = 0.0;
             for (long int j = 0; j < o; j++) {
                 dum += tmp2[q*o*v+j*v+a] * newFock[i*o+j];
             }
             Qov[q*o*v+i*v+a] = dum;
         }
     }
  }
  free(tmp2);

  // stick nQ in process environment so ccsd can know it
  Process::environment.globals["DFCC NAUX"] = (double)nQ;

  if ( memory < 8L*(o*o*v*v+o*v*nQ) ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  // allocate memory for a couple of buffers
  double*amps2 = (double*)malloc(o*o*v*v*sizeof(double));

  // build (ia|jb) integrals
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,amps2,o*v);
  free(Qov);

  if ( memory < 16L*o*o*v*v ) {
      throw PsiException("not enough memory (fno)",__FILE__,__LINE__);
  }

  double*amps1 = (double*)malloc(o*o*v*v*sizeof(double));

  // build mp2 amplitudes for mp2 density
  long int ijab = 0;
  emp2 = 0.0;
  double emp2_os = 0.0;
  double emp2_ss = 0.0;
  for (long int a=o; a<o+v; a++){
      double da = F[a];
      for (long int b=o; b<o+v; b++){
          double dab = da + F[b];
          for (long int i=0; i<o; i++){
              double dabi = dab - F[i];
              for (long int j=0; j<o; j++){
                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  double dijab = dabi-F[j];
                  amps1[ijab++] = - amps2[iajb]/dijab;
                  emp2_os -= amps2[iajb] * amps2[iajb]/dijab;
                  emp2_ss -= amps2[iajb] * (amps2[iajb] - amps2[j*o*v*v+(a-o)*o*v+i*v+(b-o)])/dijab;
              }
          }
      }
  }
  emp2 = emp2_os + emp2_ss;

  fprintf(outfile,"        Doubles contribution to MP2 energy in full space: %20.12lf\n",emp2);
  double emp2_s = 0.0;
  for (long int a = 0; a < v; a++) {
      double da = F[a+ndoccact];
      for (long int i = 0; i < o; i++) {
          double dia = da - F[i];
          double dum = - dfFock[(i+nfzc)*nmo+(a+ndocc)]/dia;
          emp2_s -= 2.0 * dum * dfFock[(i+nfzc)*nmo+(a+ndocc)];
      }
  }
  fprintf(outfile,"        Singles contribution to MP2 energy in full space: %20.12lf\n",emp2_s);
  fprintf(outfile,"\n");

  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + Process::environment.globals["SCF TOTAL ENERGY"];

  ijab = 0;
  for (long int a=o; a<o+v; a++){
      for (long int b=o; b<o+v; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){
                  long int ijba = (b-o)*o*o*v+(a-o)*o*o+i*o+j;
                  amps2[ijab] = 2.0*amps1[ijab] - amps1[ijba];
                  ijab++;
              }
          }
      }
  }

  // build ab block of the density:
  F_DGEMM('t','n',v,v,v*o*o,2.0,amps1,v*o*o,amps2,v*o*o,0.0,Dab,v);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);

  // reorder transformation matrix:
  for (long int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("OCC_TOLERANCE");
  nvirt_no = 0;
  for (long int i=0; i<v; i++) if (eigvalDab[i]>cutoff) nvirt_no++;

  fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
  fprintf(outfile,"\n");
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5li\n",v);
  fprintf(outfile,"        Number of virtual orbitals in truncated space: %5li\n",nvirt_no);
  fprintf(outfile,"\n");

  // transform Fock matrix to MP2 NO basis
  memset((void*)newFock,'\0',v*v*sizeof(double));
  F_DCOPY(v,F+ndoccact,1,newFock,v+1);
  F_DGEMM('n','n',v,nvirt_no,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',nvirt_no,nvirt_no,v,1.0,temp,v,Dab,v,0.0,newFock,nvirt_no);

  // diagonalize new Fock matrix for semi-canonical orbitals
  Diagonalize(nvirt_no,newFock,neweps);

  // construct full mo -> no transformation matrix
  F_DGEMM('n','n',v,nvirt_no,nvirt_no,1.0,temp,v,newFock,nvirt_no,0.0,Dab,v);

  // put orbital energies back in F - doesn't matter in this implementation
  F_DCOPY(nvirt_no,neweps,1,F+ndoccact,1);

  // free memory before using libtrans
  free(temp);
  free(neweps);
  free(amps1);
  free(amps2);
  free(eigvalDab);
  free(newFock);

  // change number of frozen virtual orbitals
  frzvpi_[0] = nfzv+(nvirt-nvirt_no);

  // modify c matrix
  ModifyCa(Dab);

  free(Dab);
}

void DFFrozenNO::ModifyCa(double*Dab){

  long int v = nvirt;

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();

  boost::shared_ptr<Matrix> Caomo = ref->Ca();

  double**Capointer = Caomo->pointer();

  // so->no
  double*temp = (double*)malloc(nso*nvirt_no*sizeof(double));
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<nvirt_no; j++){
          double dum = 0.0;
          for (long int k=0; k<v; k++){
              dum += Capointer[i][ndocc+k] * Dab[j*v+k];
          }
          temp[i*nvirt_no+j] = dum;
      }
  }
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<nvirt_no; j++){
          Capointer[i][ndocc+j] = temp[i*nvirt_no+j];
      }
  }
  free(temp);
}
void DFFrozenNO::ModifyCa_occ(double*Dij){

  long int o = ndoccact;

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();

  boost::shared_ptr<Matrix> Caomo = ref->Ca();

  double**Capointer = Caomo->pointer();

  // so->no
  double*temp = (double*)malloc(nso*o*sizeof(double));
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<o; j++){
          double dum = 0.0;
          for (long int k=0; k<o; k++){
              dum += Capointer[i][nfzc+k] * Dij[j*o+k];
          }
          temp[i*o+j] = dum;
      }
  }
  for (long int i=0; i<nso; i++){
      for (long int j=0; j<o; j++){
          Capointer[i][nfzc+j] = temp[i*o+j];
      }
  }
  free(temp);
}

void DFFrozenNO::BuildFock(long int nQ,double*Qso,double*F) {

    double ** Cap = Ca()->pointer();

    // Transform Qso to MO basis:
    double * tmp = (double*)malloc(nso*nso*nQ*sizeof(double));
    F_DCOPY(nso*nso*nQ,Qso,1,tmp,1);
    F_DGEMM('n','n',nmo,nso*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        for (long int mu = 0; mu < nso; mu++) {
            F_DCOPY(nmo,Qso+q*nso*nmo+mu*nmo,1,tmp+q*nso*nmo+mu,nmo);
        }
    }
    F_DGEMM('n','n',nmo,nmo*nQ,nso,1.0,&Cap[0][0],nmo,tmp,nso,0.0,Qso,nmo);

    // build Fock matrix:

    // transform H
    // one-electron integrals
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix H = mints->so_kinetic();
    H->add(mints->so_potential());

    long int max = nQ > nso*nso ? nQ : nso*nso;
    double * temp2 = (double*)malloc(max*sizeof(double));
    double * temp3 = (double*)malloc(nso*nso*sizeof(double));
    double * h     = (double*)malloc(nmo*nmo*sizeof(double));
    double ** hp   = H->pointer();

    F_DGEMM('n','t',nso,nmo,nso,1.0,&hp[0][0],nso,&Cap[0][0],nmo,0.0,temp2,nso);
    F_DGEMM('n','n',nmo,nmo,nso,1.0,&Cap[0][0],nmo,temp2,nso,0.0,h,nmo);

    // build Fock matrix:  sum_k (Q|kk)
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        double dum = 0.0;
        for (long int k = 0; k < ndocc; k++) {
            dum += Qso[q*nmo*nmo + k*nmo + k];
        }
        temp2[q] = 2.0 * dum;
    }

    // use temp and Qso as storage for Qmo( q, r, k)
    #pragma omp parallel for schedule (static)
    for (long int q = 0; q < nQ; q++) {
        for (long int r = 0; r < nmo; r++) {
            for (long int k = 0; k < ndocc; k++) {
                tmp[q*nmo*ndocc+k*nmo+r]  = Qso[q*nmo*nmo+k*nmo+r];
            }
        }
    }
    // I(r,s) = sum q k (q|ks)(q|rk)
    F_DGEMM('n','t',nmo,nmo,nQ*ndocc,1.0,tmp,nmo,tmp,nmo,0.0,temp3,nmo);

    // Fock matrix
    #pragma omp parallel for schedule (static)
    for (long int i = 0; i < nmo; i++) {
        for (long int j = 0; j < nmo; j++) {
            double dum = h[i*nmo+j] - temp3[i*nmo+j];
            dum += F_DDOT(nQ,temp2,1,Qso + i*nmo + j , nmo*nmo);
            F[i*nmo+j] = dum;
        }
    }

    free(h);
    free(tmp);
    free(temp2);
    free(temp3);
}

}} // end of namespaces
