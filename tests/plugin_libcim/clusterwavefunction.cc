#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>

using namespace psi;
using namespace boost;

namespace psi{

/*
 * change reference wavefunction to pass to correlated calculation.  
 * all that changes is the number of frozen occupied and virtual orbitals.
 */
void CIM::ClusterWavefunction(int cluster){
  int nfzc_  = nfzc;
  int nfzv_  = nfzv;

  int nactiveo   = ncentral[cluster]+nmodomain[cluster]+nenv[cluster];
  int ninactiveo = ndocc - nfzc - nactiveo;
  int nactivev   = nvirt_;
  int ninactivev = nvirt-nvirt_;

  nfzc_  += ninactiveo;
  nfzv_  += ninactivev;

  frzcpi_[0] = nfzc_;
  frzvpi_[0] = nfzv_;
printf("hey %5i %5i %5i %5i\n",nactiveo,ninactiveo,frzcpi()[0],frzcpi()[0]);
}

}
