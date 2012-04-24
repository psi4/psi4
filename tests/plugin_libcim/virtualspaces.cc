#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"cim.h"
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<lib3index/dftensor.h>
#include<lib3index/3index.h>

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
  //MP2Density(Dab,cluster);
  OppositeSpinMP2Density(Dab,cluster);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);

  // reorder transformation matrix (so eigenvals are large->small)
  double*temp    = (double*)malloc(v*v*sizeof(double));
  for (int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options_.get_double("OCC_TOLERANCE");
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
  double**Fock = boys->Fock->pointer();

  // eigenvalues (actually just the diagonals of F)
  shared_ptr<Vector> evals_aocc_ (new Vector("Epsilon (Active Occupied)",o));
  shared_ptr<Vector> evals_avir_ (new Vector("Epsilon (Active Virtual)",v));
  double *evals_aocc_pointer = evals_aocc_->pointer();
  double *evals_avir_pointer = evals_avir_->pointer();
  for (int i=0; i<o; i++) evals_aocc_pointer[i] = Fock[i+nfzc][i+nfzc];
  for (int a=0; a<v; a++) evals_avir_pointer[a] = Fock[a+ndocc][a+ndocc];

  // laplace denominator:
  boost::shared_ptr<Denominator> 
      denom(Denominator::buildDenominator("LAPLACE", evals_aocc_, evals_avir_,options_.get_double("DENOMINATOR_DELTA")));

  int nW = denom->nvector();
  double *Qw1 = (double*)malloc(o*v*nQ*sizeof(double));
  double *Qw2 = (double*)malloc(o*v*nQ*sizeof(double));
  SharedMatrix tau = denom->denominator();
  double **taup = tau->pointer();

  double *I = (double*)malloc(nQ*nQ*sizeof(double));
  double *Ip = (double*)malloc(o*v*nQ*sizeof(double));

  memset((void*)Dab,'\0',v*v*sizeof(double));
  for (int w1=0; w1<nW; w1++){
      F_DCOPY(o*v*nQ,Qov[0],1,Qw1,1);
      for (int q=0; q<nQ; q++){
          for (int ia=0; ia<o*v; ia++){
              Qw1[q*o*v+ia] *= taup[w1][ia];
          }
      }
      for (int w2=0; w2<nW; w2++){
          F_DCOPY(o*v*nQ,Qov[0],1,Qw2,1);
          for (int q=0; q<nQ; q++){
              for (int ia=0; ia<o*v; ia++){
                  Qw2[q*o*v+ia] *= taup[w2][ia];
              }
          }

          // I(Q,P) = (Q|jc)(P|jc)
          F_DGEMM('t','n',nQ,nQ,o*v,1.0,Qw2,o*v,Qw1,o*v,0.0,I,nQ);
          // I'(P,ia) = I(Q,P)(ia|Q)
          F_DGEMM('n','t',o*v,nQ,nQ,1.0,Qw1,o*v,I,nQ,0.0,Ip,o*v);
          // I''(iab) = I'(P,ia)(ib|P)
          for (int a=0; a<v; a++){
              for (int b=a; b<v; b++){
                  double dum = 0.0;
                  for (int i=0; i<o; i++){
                      if (domain[cluster][i]==isempty) continue;
                      dum += F_DDOT(nQ,Ip+i*v+a,o*v,Qw2+i*v+b,o*v);
                  }
                  Dab[a*v+b] += dum;
                  if (a!=b) Dab[b*v+a] += dum;
              }
          }
      }
  }
  for (int a=0; a<v; a++){
      for (int b=0; b<v; b++){
          Dab[a*v+b] *= 2.0*1.3;
      }
  }
  free(I);
  free(Ip);
  free(Qw1);
  free(Qw2);

  /*double*t2 = (double*)malloc(o*o*v*v*sizeof(double));
  memset((void*)t2,'\n',o*o*v*v*sizeof(double));
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
  free(t2);*/
}



}
