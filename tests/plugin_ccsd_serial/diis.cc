#include"psi4-dec.h"
#include<psifiles.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>

#include"ccsd.h"
#include"sort.h"
#include"blas.h"

using namespace psi;

/*================================================================

   diis functions

================================================================*/

namespace psi{

void CoupledCluster::DIIS(double*c,long int nvec,long int n){
  long int i,j,k;
  doublereal sum,dum;
  integer*ipiv,nvar;
  nvar = nvec+1;
  doublereal*A = (doublereal*)malloc(sizeof(doublereal)*nvar*nvar);
  doublereal*B = (doublereal*)malloc(sizeof(doublereal)*nvar);
  memset((void*)A,'\0',nvar*nvar*sizeof(double));
  memset((void*)B,'\0',nvar*sizeof(double));
  B[nvec] = -1.;
  ipiv = (integer*)malloc(nvar*sizeof(integer));

  char*evector=(char*)malloc(1000*sizeof(char));

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_EVEC,PSIO_OPEN_OLD);

  for (i=0; i<nvec; i++){
      sprintf(evector,"evector%li",i+1);
      psio->read_entry(PSIF_EVEC,evector,(char*)&tempt[0],n*sizeof(double));
      for (j=i+1; j<nvec; j++){
          sprintf(evector,"evector%li",j+1);
          psio->read_entry(PSIF_EVEC,evector,(char*)&tempv[0],n*sizeof(double));
          sum = F_DDOT(n,tempt,1,tempv,1);
          A[j*nvar+i] = sum;
          A[i*nvar+j] = sum;
      }
      A[i*nvar+i] = F_DDOT(n,tempt,1,tempt,1);
  }
  j = nvec;
  for (i=0; i<nvar; i++){
      A[j*nvar+i] = -1.0;
      A[i*nvar+j] = -1.0;
  }
  A[nvar*nvar-1] = 0.;
  psio->close(PSIF_EVEC,1);
  free(evector);

  integer nrhs,lda,ldb,info;
  nrhs = 1;
  lda = ldb = nvar;
  info = 0;
  DGESV(nvar,nrhs,A,lda,ipiv,B,ldb,info);
  F_DCOPY(nvec,B,1,c,1);

  free(A);
  free(B);
  free(ipiv);
  psio.reset();
}
void CoupledCluster::DIISOldVector(long int iter,int diis_iter,int replace_diis_iter){
  long int j,o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector=(char*)malloc(1000*sizeof(char));

  if (diis_iter<=maxdiis && iter<=maxdiis){
     sprintf(oldvector,"oldvector%i",diis_iter);
  }
  else{
     sprintf(oldvector,"oldvector%i",replace_diis_iter);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0)
     psio->open(PSIF_OVEC,PSIO_OPEN_NEW);
  else
     psio->open(PSIF_OVEC,PSIO_OPEN_OLD);

  psio_address addr;
  addr = PSIO_ZERO;

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = integrals;
  }

  psio->write(PSIF_OVEC,oldvector,(char*)&tb[0],arraysize*sizeof(double),addr,&addr);
  psio->write(PSIF_OVEC,oldvector,(char*)&t1[0],o*v*sizeof(double),addr,&addr);
  psio->close(PSIF_OVEC,1);
  psio.reset();

  free(oldvector);
}
double CoupledCluster::DIISErrorVector(int diis_iter,int replace_diis_iter,int iter){
  double nrm;
  long int i,j,o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*evector   = (char*)malloc(1000*sizeof(char));
  if (diis_iter<=maxdiis && iter<=maxdiis){
     sprintf(evector,"evector%i",diis_iter);
  }
  else{
     sprintf(evector,"evector%i",replace_diis_iter);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0)
     psio->open(PSIF_EVEC,PSIO_OPEN_NEW);
  else
     psio->open(PSIF_EVEC,PSIO_OPEN_OLD);

  nrm = F_DNRM2(arraysize+o*v,tempv,1);
  psio->write_entry(PSIF_EVEC,evector,(char*)&tempv[0],(arraysize+o*v)*sizeof(double));

  psio->close(PSIF_EVEC,1);
  psio.reset();

  free(evector);

  // return convergence
  return nrm;
}
void CoupledCluster::DIISNewAmplitudes(int diis_iter){
  long int o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector;
  oldvector=(char*)malloc(1000*sizeof(char));

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_OVEC,PSIO_OPEN_OLD);

  psio_address addr;

  if (t2_on_disk){
     tb = integrals;
  }

  memset((void*)tb,'\0',arraysize*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));

  int max = diis_iter;
  if (max > maxdiis) max = maxdiis;

  for (long int j=1; j<=max; j++){
      addr = PSIO_ZERO;
      sprintf(oldvector,"oldvector%li",j);
      psio->read(PSIF_OVEC,oldvector,(char*)&tempt[0],arraysize*sizeof(double),addr,&addr);
      F_DAXPY(arraysize,diisvec[j-1],tempt,1,tb,1);
      psio->read(PSIF_OVEC,oldvector,(char*)&tempt[0],o*v*sizeof(double),addr,&addr);
      F_DAXPY(o*v,diisvec[j-1],tempt,1,t1,1);
  }
  psio->close(PSIF_OVEC,1);
  free(oldvector);

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_T2,"t2",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }

  psio.reset();
}

}
