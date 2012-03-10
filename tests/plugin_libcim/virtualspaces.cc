#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<lib3index/dftensor.h>

using namespace psi;
using namespace boost;

namespace psi{

/*
 * build virtual spaces based on MP2 OPDM
 */
void CIM::VirtualSpaces(int cluster,int clusternum){

  fprintf(outfile,"  ==> Build Virtual Space for Cluster {%i} <==\n",clusternum);
  fprintf(outfile,"\n");

  int o = ndoccact;
  int v = nvirt;

  // build mp2 density using canonical formulas
  double*Dab=(double*)malloc(v*v*sizeof(double));
  MP2Density(Dab,cluster);
  //OppositeSpinMP2Density(Dab,cluster);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);

  // reorder transformation matrix (so eigenvals are large->small)
  double*temp    = (double*)malloc(v*v*sizeof(double));
  for (int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("THRESH3");
  nvirt_ = 0;
  for (int i=0; i<v; i++) if (eigvalDab[i]>=cutoff) nvirt_++;

  fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
  fprintf(outfile,"\n");
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5i\n",v);
  fprintf(outfile,"        Number of virtual orbitals in truncated space: %5i\n",nvirt_);
  fprintf(outfile,"\n");

  // transform Fock matrix to NO basis ( for canonicalization in QuasiCanonicalOrbitals() )
  double**Fock = boys->Fock->pointer();
  double*newFock = (double*)malloc(nso*nmo*sizeof(double));
  memset((void*)newFock,'\0',v*v*sizeof(double));
  for (int a=0; a<v; a++){
      for (int b=0; b<v; b++){
          newFock[a*v+b] = Fock[a+ndocc][b+ndocc];
      }
  }
  F_DGEMM('n','n',v,v,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',v,v,v,1.0,temp,v,Dab,v,0.0,newFock,v);

  // overwrite localFock virtual-virtual block
  Fock = localFock->pointer();
  for (int a=0; a<nvirt_; a++){
      for (int b=0; b<nvirt_; b++){
          Fock[a+ndocc][b+ndocc] = newFock[a*v+b];
      }
  }

  // modify C matrix to map SO->NO
  Fock = boys->Clmo->pointer();
  for (int i=0; i<nso; i++){
      for (int a=0; a<v; a++){
          double dum = 0.0;
          for (int b=0; b<v; b++){
              dum += temp[a*v+b] * Fock[i][ndocc+b];
          }
          newFock[i*nmo+ndocc+a] = dum;
      }
  }
  // overwrite localClmo virtual-virtual block
  Fock = localClmo->pointer();
  for (int i=0; i<nso; i++){
      for (int a=0; a<v; a++){
          Fock[i][ndocc+a] = newFock[i*nmo+ndocc+a];
      }
  }
  
  free(temp);
  free(Dab);
  free(newFock);
}

void CIM::MP2Density(double*Dab,int cluster){
  int o = ndoccact;
  int v = nvirt;
  double*t2 = (double*)malloc(o*o*v*v*sizeof(double));
  memset((void*)t2,'\n',o*o*v*v*sizeof(double));
  double**Fock = boys->Fock->pointer();

  F_DGEMM('n','t',o*v,o*v,nQ,1.0,&Qov[0][0],o*v,&Qov[0][0],o*v,0.0,t2,o*v);
  for (int a=0; a<v; a++){
      double da = -Fock[a+ndocc][a+ndocc];
      for (int b=0; b<v; b++){
          double dab = da-Fock[b+ndocc][b+ndocc];
          for (int i=0; i<o; i++){
              double dabi = dab+Fock[i+nfzc][i+nfzc];
              for (int j=0; j<o; j++){
                  double dabij = dabi+Fock[j+nfzc][j+nfzc];
                  t2[i*o*v*v+a*o*v+j*v+b] /= dabij;
              }
          }
      }
  }

  // build Dab({I})
  for (int a=0; a<v; a++){
      for (int b=a; b<v; b++){
          double dum = 0.0;
          for (int i=0; i<o; i++){
              if (domain[cluster][i]==isempty) continue;
              for (int j=0; j<o; j++){
                  for (int c=0; c<v; c++){
                      dum += t2[i*o*v*v+a*o*v+j*v+c] * (2.0*t2[i*o*v*v+b*o*v+j*v+c] - t2[i*o*v*v+c*o*v+j*v+b]);
                  }
              }
          }
          Dab[a*v+b] = 2.0 * dum;
          Dab[b*v+a] = 2.0 * dum;
      }
  }
  free(t2);
}

void CIM::OppositeSpinMP2Density(double*Dab,int cluster){
  int o = ndoccact;
  int v = nvirt;
  double*t2 = (double*)malloc(o*o*v*v*sizeof(double));
  memset((void*)t2,'\n',o*o*v*v*sizeof(double));
  double**Fock = boys->Fock->pointer();

  F_DGEMM('n','t',o*v,o*v,nQ,1.0,&Qov[0][0],o*v,&Qov[0][0],o*v,0.0,t2,o*v);
  for (int a=0; a<v; a++){
      double da = -Fock[a+ndocc][a+ndocc];
      for (int b=0; b<v; b++){
          double dab = da-Fock[b+ndocc][b+ndocc];
          for (int i=0; i<o; i++){
              double dabi = dab+Fock[i+nfzc][i+nfzc];
              for (int j=0; j<o; j++){
                  double dabij = dabi+Fock[j+nfzc][j+nfzc];
                  t2[i*o*v*v+a*o*v+j*v+b] /= dabij;
              }
          }
      }
  }

  // build Dab({I})
  for (int a=0; a<v; a++){
      for (int b=a; b<v; b++){
          double dum = 0.0;
          for (int i=0; i<o; i++){
              if (domain[cluster][i]==isempty) continue;
              for (int j=0; j<o; j++){
                  for (int c=0; c<v; c++){
                      dum += t2[i*o*v*v+a*o*v+j*v+c] * t2[i*o*v*v+b*o*v+j*v+c];
                  }
              }
          }
          Dab[a*v+b] = 2.0 * dum;
          Dab[b*v+a] = 2.0 * dum;
      }
  }
  free(t2);
}

}
