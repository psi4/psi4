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
#include<../src/bin/cepa/blas.h>
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>

using namespace psi;
using namespace cepa;
using namespace boost;
namespace psi{namespace qci{
void SortOVOV(struct iwlbuf *Buf,int nfzc,int nfzv,int norbs,int ndoccact,int nvirt);
}}

namespace psi{namespace qci{

FrozenNO::FrozenNO(Options&options):
  Wavefunction(Process::environment.options, _default_psio_lib_),options_(options)
{
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.wavefunction();
  copy(ref);

  if (ref.get() == NULL){
     throw PsiException("no wavefunction?",__FILE__,__LINE__);
  }
  nirreps = ref->nirrep();
  if (nirreps>1){
     throw PsiException("FrozenNO requires symmetry c1 (for now!)",__FILE__,__LINE__);
  }
  int * sorbs = ref->nsopi();
  int * orbs  = ref->nmopi();
  int * docc  = ref->doccpi();
  int * fzc   = ref->frzcpi();
  int * fzv   = ref->frzvpi();
  nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
  for (long int h=0; h<nirreps; h++){
      nfzc   += fzc[h];
      nfzv   += fzv[h];
      nso    += sorbs[h];
      nmo    += orbs[h];
      ndocc  += docc[h];
  }
  ndoccact = ndocc - nfzc;
  nvirt  = nmo - ndocc;

  // NO time
  NaturalOrbitals();
}
FrozenNO::~FrozenNO()
{
}
// use this function to return the mp2 energy in the full basis.
double FrozenNO::compute_energy(){
  return emp2;
}

/*
 * build natural orbitals and transform TEIs
 */
void FrozenNO::NaturalOrbitals(){
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

  // allocate memory for a couple of buffers (TODO: out-of-core)
  long int o = ndoccact;
  long int v = nvirt;
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
                  amps2[ijab] = 2.0*amps1[ijab++] - amps1[ijba];
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

}} // end of namespaces
