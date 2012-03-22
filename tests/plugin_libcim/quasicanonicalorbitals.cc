#include<psifiles.h>
#include<libpsio/psio.hpp>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>


using namespace psi;
using namespace boost;

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

/* 
 * define quasicanonical orbitals for cluster {i}
 */
void CIM::QuasiCanonicalOrbitals(int cluster){

  double **C = boys->Clmo->pointer();
  double **F = boys->Fock->pointer();

  // how many times does each central orbital occur?
  int*ntimes = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++) ntimes[i] = 0;
  // which orbtals are central in this cluster?
  for (int j=0; j<ndoccact; j++){
      if (central[cluster][j]!=isempty) ntimes[j]=1;
  }
  // which orbtals are central in all other clusters?
  for (int j=0; j<maxndomains; j++){
      if (skip[j]) continue;
      if (cluster==j) continue;
      for (int k=0; k<ndoccact; k++){
          if (central[j][k]==isempty) continue;
          if (ntimes[k]!=0) ntimes[k]++;
      }
  }
  // reorder ntimes
  int nact = 0;
  int*tempntimes = (int*)malloc(ndoccact*sizeof(int));
  for (int i=0; i<ndoccact; i++) tempntimes[i] = ntimes[i];
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]!=isempty){
         tempntimes[nact++] = ntimes[i];
      }
  }
  for (int i=0; i<ndoccact; i++) ntimes[i] = tempntimes[i];
  free(tempntimes);

  // map mo->lmo transformation and lmo Fock matrix s.t. inactive orbitals come first
  int ninact = 0;
  nact = 0;

  double**localC_pointer = localClmo->pointer();
  double**localF_pointer = localFock->pointer();

  // reorder C
  for (int i=0; i<nfzc; i++){
      for (int j=0; j<nso; j++){
          localC_pointer[j][i] = C[j][i];
      }
  }
  ninact = nfzc;
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]==isempty){
         for (int j=0; j<nso; j++){
             localC_pointer[j][ninact] = C[j][i+nfzc];
         }
         ninact++;
      }
  }
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]!=isempty){
         for (int j=0; j<nso; j++){
             localC_pointer[j][ninact+nact] = C[j][i+nfzc];
         }
         nact++;
      }
  }

  // copy all nvirt, not just nvirt_.
  //for (int i=0; i<nvirt; i++){
  //    for (int j=0; j<nso; j++){
  //        localC_pointer[j][i+ninact+nact] = C[j][i+ndocc];
  //    }
  //}

  // reorder F (1st index)
  double*tempF = (double*)malloc(nmo*nmo*sizeof(double));
  nact = 0;
  for (int i=0; i<nfzc; i++){
      for (int j=0; j<nmo; j++){
          tempF[i*nmo+j] = F[i][j];
      }
  }
  ninact = nfzc;
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]==isempty){
         for (int j=0; j<nmo; j++){
             tempF[ninact*nmo+j] = F[i+nfzc][j];
         }
         ninact++;
      }
  }
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]!=isempty){
         for (int j=0; j<nmo; j++){
             tempF[(ninact+nact)*nmo+j] = F[i+nfzc][j];
         }
         nact++;
      }
  }
  // don't copy Fvirt from F, localFock was filled in VirtualOrbitals()
  // TODO this won't be nvirt once i build the virtual space for each cluster
  //for (int i=0; i<nvirt; i++){
  //    for (int j=0; j<nmo; j++){
  //        tempF[(i+ninact+nact)*nmo+j] = F[i+ndocc][j];
  //    }
  //}


  // reorder F (2nd index)
  nact = 0;
  for (int i=0; i<nfzc; i++){
      for (int j=0; j<nmo; j++){
          localF_pointer[j][i] = tempF[j*nmo+i];
      }
  }
  ninact = nfzc;
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]==isempty){
         for (int j=0; j<nmo; j++){
             localF_pointer[j][ninact] = tempF[j*nmo+i+nfzc];
         }
         ninact++;
      }
  }
  for (int i=0; i<ndoccact; i++){
      if (domain[cluster][i]!=isempty){
         for (int j=0; j<nmo; j++){
             localF_pointer[j][ninact+nact] = tempF[j*nmo+i+nfzc];
         }
         nact++;
      }
  }
  // copy all nvirt, not just nvirt_.
  //for (int i=0; i<nvirt; i++){
  //    for (int j=0; j<nmo; j++){
  //        localF_pointer[j][i+ninact+nact] = tempF[j*nmo+i+ndocc];
  //    }
  //}
  free(tempF);

  // diagonalize active occupied block of Fock matrix
  double*eigval=(double*)malloc(nmo*sizeof(double));
  double*mat=(double*)malloc(nmo*nmo*sizeof(double));
  double*temp1=(double*)malloc(nso*nmo*sizeof(double));
  for (int i=0; i<nact; i++){
      for (int j=0; j<nact; j++){
          mat[i*nact+j] = localF_pointer[ninact+i][ninact+j];
          //mat[i*nact+j] = F[ninact+i][ninact+j];
      }
  }
  Diagonalize(nact,mat,eigval);
  double**Rii_pointer = Rii->pointer();
  for (int i=0; i<nact; i++){
      for (int j=0; j<nact; j++){
          Rii_pointer[i][j] = mat[i*nact+j];
      }
  }

  // replace eps from wavefunction:
  double*eps = epsilon_a()->pointer();
  // save a copy of eps to re-replace after CC runs. this is a stupid hack.
  // apparently copy(ref) doesn't clone the reference, it just points me 
  // toward it.  so when i modify epsilon_a() in cim, i modify it in 
  // the reference
  epsSave = (double*)malloc(nmo*sizeof(double));
  F_DCOPY(nmo,eps,1,epsSave,1);

  for (int i=ninact; i<ndocc; i++) eps[i] = eigval[i-ninact];

  // transform so->lmo to so->quasicanconical lmo
  for (int i=0; i<nso; i++){
      for (int j=0; j<nact; j++){
          double dum = 0.0;
          for (int k=0; k<nact; k++){
              dum += localC_pointer[i][ninact+k] * mat[j*nact+k];
          }
          temp1[i*nact+j] = dum;
      }
  }
  for (int i=0; i<nso; i++){
      for (int j=0; j<nact; j++){
          localC_pointer[i][ninact+j] = temp1[i*nact+j];
      }
  }

  // diagonalize active virtual block of Fock matrix
  for (int i=0; i<nvirt_; i++){
      for (int j=0; j<nvirt_; j++){
          mat[i*nvirt_+j] = localF_pointer[ndocc+i][ndocc+j];
      }
  }
  Diagonalize(nvirt_,mat,eigval);

  // replace eps from wavefunction:
  for (int i=ndocc; i<ndocc+nvirt_; i++) eps[i] = eigval[i-ndocc];
 
  // transform so->lmo to so->quasicanconical lmo
  for (int i=0; i<nso; i++){
      for (int j=0; j<nvirt_; j++){
          double dum = 0.0;
          for (int k=0; k<nvirt_; k++){
              dum += localC_pointer[i][ndocc+k] * mat[j*nvirt_+k];
          }
          temp1[i*nvirt_+j] = dum;
      }
  }
  for (int i=0; i<nso; i++){
      for (int j=0; j<nvirt_; j++){
          localC_pointer[i][ndocc+j] = temp1[i*nvirt_+j];
      }
  }
  free(temp1);

  // factor for scaling contribution from each central orbital
  for (int i=0; i<nact; i++){
      if (ntimes[i]==0) centralfac[i] = 0.0;
      else              centralfac[i] = 1.0/ntimes[i];
  }
  free(ntimes);
}

}
