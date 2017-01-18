/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include"ccsd.h"
#include"blas.h"
#include"psi4/libqt/qt.h"

using namespace psi;

/*================================================================

   diis functions

================================================================*/

namespace psi{ namespace fnocc{

void CoupledCluster::DIIS(double*c,long int nvec,long int n,int replace_diis_iter){
  integer nvar      = nvec+1;
  integer    * ipiv = (integer*)malloc(nvar*sizeof(integer));
  doublereal * temp = (doublereal*)malloc(sizeof(doublereal)*maxdiis*maxdiis);
  doublereal * A    = (doublereal*)malloc(sizeof(doublereal)*nvar*nvar);
  doublereal * B    = (doublereal*)malloc(sizeof(doublereal)*nvar);
  memset((void*)A,'\0',nvar*nvar*sizeof(double));
  memset((void*)B,'\0',nvar*sizeof(double));
  B[nvec] = -1.;

  char*evector=(char*)malloc(1000*sizeof(char));

  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_EVEC,PSIO_OPEN_OLD);

  // add row to matrix, don't build the whole thing.
  psio->read_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis*maxdiis*sizeof(double));
  for (long int i = 0; i < nvec; i++){
      for (long int j = 0; j < nvec; j++){
          A[i*nvar+j] = temp[i*maxdiis+j];
      }
  }

  if (nvec <= 3) {
      for (long int i = 0; i < nvec; i++) {
          sprintf(evector,"evector%li",i+1);
          psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&tempt[0],n*sizeof(double));
          for (long int j = i; j < nvec; j++){
              sprintf(evector,"evector%li",j+1);
              psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&tempv[0],n*sizeof(double));
              double sum  = C_DDOT(n,tempt,1,tempv,1);
              A[i*nvar+j] = sum;
              A[j*nvar+i] = sum;
          }
      }
  }else {
      long int i;
      if (nvec<=maxdiis && iter<=maxdiis){
          i = nvec - 1;
      }
      else{
          i = replace_diis_iter - 1;
      }
      sprintf(evector,"evector%li",i+1);
      psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&tempt[0],n*sizeof(double));
      for (long int j = 0; j < nvec; j++){
          sprintf(evector,"evector%li",j+1);
          psio->read_entry(PSIF_DCC_EVEC,evector,(char*)&tempv[0],n*sizeof(double));
          double sum  = C_DDOT(n,tempt,1,tempv,1);
          A[i*nvar+j] = sum;
          A[j*nvar+i] = sum;
      }
  }

  long int j = nvec;
  for (long int i = 0; i < nvar; i++){
      A[j*nvar+i] = -1.0;
      A[i*nvar+j] = -1.0;
  }
  A[nvar*nvar-1] = 0.;

  // save matrix for next iteration
  for (long int i = 0; i < nvec; i++){
      for (long int j = 0; j < nvec; j++){
          temp[i*maxdiis+j] = A[i*nvar+j];
      }
  }
  psio->write_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis*maxdiis*sizeof(double));
  free(temp);
  psio->close(PSIF_DCC_EVEC,1);
  free(evector);

  integer nrhs,lda,ldb,info;
  nrhs = 1;
  lda = ldb = nvar;
  info = 0;
  DGESV(nvar,nrhs,A,lda,ipiv,B,ldb,info);
  C_DCOPY(nvec,B,1,c,1);

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

  std::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0) {
     psio->open(PSIF_DCC_OVEC,PSIO_OPEN_NEW);
  }else {
     psio->open(PSIF_DCC_OVEC,PSIO_OPEN_OLD);
  }

  psio_address addr;
  addr = PSIO_ZERO;

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }

  psio->write(PSIF_DCC_OVEC,oldvector,(char*)&tb[0],arraysize*sizeof(double),addr,&addr);
  psio->write(PSIF_DCC_OVEC,oldvector,(char*)&t1[0],o*v*sizeof(double),addr,&addr);
  psio->close(PSIF_DCC_OVEC,1);
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

  std::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0) {
     psio->open(PSIF_DCC_EVEC,PSIO_OPEN_NEW);
     double * temp = (double*)malloc(maxdiis*maxdiis*sizeof(double));
     memset((void*)temp,'\0',maxdiis*maxdiis*sizeof(double));
     psio->write_entry(PSIF_DCC_EVEC,"error matrix",(char*)&temp[0],maxdiis*maxdiis*sizeof(double));
     free(temp);
  }
  else {
     psio->open(PSIF_DCC_EVEC,PSIO_OPEN_OLD);
  }

  nrm = C_DNRM2(arraysize+o*v,tempv,1);
  psio->write_entry(PSIF_DCC_EVEC,evector,(char*)&tempv[0],(arraysize+o*v)*sizeof(double));

  psio->close(PSIF_DCC_EVEC,1);
  psio.reset();

  free(evector);

  // return convergence
  return nrm;
}
void CoupledCluster::DIISNewAmplitudes(int diis_iter,int&replace_diis_iter){
  long int o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector;
  oldvector=(char*)malloc(1000*sizeof(char));

  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_OVEC,PSIO_OPEN_OLD);

  psio_address addr;

  if (t2_on_disk){
     tb = integrals;
  }

  memset((void*)tb,'\0',arraysize*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));

  long int max = diis_iter;
  if (max > maxdiis) max = maxdiis;

  double min = 1.e9;
  for (long int j=1; j<=max; j++){
      addr = PSIO_ZERO;
      sprintf(oldvector,"oldvector%li",j);
      psio->read(PSIF_DCC_OVEC,oldvector,(char*)&tempt[0],arraysize*sizeof(double),addr,&addr);
      C_DAXPY(arraysize,diisvec[j-1],tempt,1,tb,1);
      psio->read(PSIF_DCC_OVEC,oldvector,(char*)&tempt[0],o*v*sizeof(double),addr,&addr);
      C_DAXPY(o*v,diisvec[j-1],tempt,1,t1,1);
      //if ( fabs( diisvec[j-1] ) < min ) {
      //    min = fabs( diisvec[j-1] );
      //    replace_diis_iter = j;
      //}
  }
  psio->close(PSIF_DCC_OVEC,1);
  free(oldvector);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }

  psio.reset();
}

}}
