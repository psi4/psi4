#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"boys.h"
#include"cim.h"
#include<libmints/mints.h>
#include<libmints/mintshelper.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>

using namespace psi;
using namespace boost;

namespace psi{

CIM::CIM(Options&options):
  Wavefunction(options, _default_psio_lib_)
{

  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                 Cluster-in-molecule                 *\n");
  fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);

  options_ = options;
  BuildClusters();
}
CIM::~CIM()
{}
double CIM::compute_energy(){
  return 0.0;
}

void CIM::BuildClusters(){
  /*
   * grab reference and parameters
   */
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  copy(ref);

  if (ref.get() == NULL){
     throw PsiException("no wavefunction?",__FILE__,__LINE__);
  }
  nirreps = ref->nirrep();
  if (nirreps>1){
     throw PsiException("cim requires symmetry c1",__FILE__,__LINE__);
  }
  int*sorbs = ref->nsopi();
  int*orbs  = ref->nmopi();
  int*docc  = ref->doccpi();
  int*fzc   = ref->frzcpi();
  int*fzv   = ref->frzvpi();
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

  /*
   * localize orbitals
   */ 
  boys  = boost::shared_ptr<Boys>(new Boys(options_));

  // print lmos
  //boys->Clmo->print();

  // print lmo fock matrix
  //boys->Fock->print();

  /*
   * grab lmo fock matrix
   */ 
  Fock = boys->Fock->pointer();

  /*
   * build occupied domains
   */
  OccupiedDomains();
}

}

