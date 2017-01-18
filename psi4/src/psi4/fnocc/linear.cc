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

#include "psi4/libmints/wavefunction.h"

#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psifiles.h"
#include"psi4/libqt/qt.h"
#ifdef _OPENMP
    #include<omp.h>
#endif

#include"blas.h"
#include"ccsd.h"

using namespace psi;

namespace psi{ namespace fnocc{

void CoupledCluster::CPU_t1_vmeai_linear(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int i,a,m,e,id,one=1;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  C_DAXPY(o*o*v*v,-2.0,integrals,1,tempv,1);

  for (i=0; i<o; i++){
      C_DCOPY(v,t1+i,o,tempt+i*v,1);
  }
  F_DGEMV('n',o*v,o*v,-1.0,tempv,o*v,tempt,1,0.0,integrals,1);
  for (a=0; a<v; a++){
      C_DAXPY(o,1.0,integrals+a,v,w1+a*o,1);
  }

  psio.reset();
}

void CoupledCluster::CPU_t1_vmeni_linear(CCTaskParams params){
  long int m,e,n,a,id;
  long int o=ndoccact;
  long int v=nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (a=0,id=0; a<v; a++){
      for (m=0; m<o; m++){
          for (n=0; n<o; n++){
              for (e=0; e<v; e++){
                  tempt[id++] = 2.*tb[e*v*o*o+a*o*o+m*o+n]-tb[a*v*o*o+e*o*o+m*o+n];
              }
          }
      }
  }
  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_DCC_IJAK,1);
  F_DGEMM('t','n',o,v,o*o*v,-1.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);
  psio.reset();
}

void CoupledCluster::CPU_t1_vmaef_linear(CCTaskParams params){
  long int m,e,i,f,a,id;
  long int o=ndoccact;
  long int v=nvirt;

  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (f=0,id=0; f<v; f++){
      for (m=0; m<o; m++){
          for (e=0; e<v; e++){
              for (i=0; i<o; i++){
                  tempt[id++] = 2.*tb[e*v*o*o+f*o*o+m*o+i]-tb[e*v*o*o+f*o*o+i*o+m];
              }
          }
      }
  }

  long int tilesize,lasttile,ntiles=1;
  long int ov2 = o*v*v;
  // tile v in chunks of o

  ntiles=1L;
  tilesize=v/1L;
  if (ntiles*tilesize<v) tilesize++;
  while(tilesize*ov2>maxelem){
     ntiles++;
     tilesize = v/ntiles;
     if (ntiles*tilesize<ov2) tilesize++;
  }
  lasttile = v - (ntiles-1L)*tilesize;

  psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;

  for (i=0; i<ntiles-1; i++){
      psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],tilesize*ov2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,tilesize,ov2,1.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  }
  i=ntiles-1;
  psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],lasttile*ov2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lasttile,ov2,1.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  psio->close(PSIF_DCC_ABCI3,1);
  psio.reset();

}

// a refactored version of I2p(ab,ci) that avoids ov^3 storage
void CoupledCluster::CPU_I2p_abci_refactored_term1_linear(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int a,b,c,i,j,id=0;
  long int ov2 = o*v*v;
  long int o2v = o*o*v;

  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_ABCI5,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;

  for (i=0; i<nov2tiles-1; i++){
      psio->read(PSIF_DCC_ABCI5,"E2abci5",(char*)&integrals[0],v*ov2tilesize*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,ov2tilesize,v,1.0,t1,o,integrals,v,0.0,tempt+i*ov2tilesize*o,o);
  }
  i=nov2tiles-1;
  psio->read(PSIF_DCC_ABCI5,"E2abci5",(char*)&integrals[0],v*lastov2tile*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lastov2tile,v,1.0,t1,o,integrals,v,0.0,tempt+i*ov2tilesize*o,o);
  psio->close(PSIF_DCC_ABCI5,1);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
  for (b=0; b<v; b++){
      C_DAXPY(o*o,1.0,tempt+b*v*o*o+a*o*o,1,tempv+a*v*o*o+b*o*o,1);
  }}
  for (a=0; a<v; a++){
  for (b=0; b<v; b++){
  for (i=0; i<o; i++){
      C_DAXPY(o,1.0,tempt+a*v*o*o+b*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
  }}}
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();
}

void CoupledCluster::UpdateT1_mp4(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int id=0;
  double tnew,dia,energy;

  if (iter<1){
     memset((void*)t1,'\0',o*v*sizeof(double));
     memset((void*)w1,'\0',o*v*sizeof(double));
  }
  else{
     for (i=0; i<o; i++){

         for (a=o; a<rs; a++){
             dia = -eps[i]+eps[a];
             tnew = - (w1[(a-o)*o+i])/dia;
             w1[(a-o)*o+i] = tnew;
         }
     }
  }
  // error vector for diis is in tempv:
  C_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  C_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  C_DCOPY(o*v,w1,1,t1,1);
}
void CoupledCluster::UpdateT2_mp4(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,di,dij,dija,energy;
  long int iajb,jaib,ijab=0;
  //double energy = 0.0;
  // we still have the residual in memory in tempv
  //psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));

  // compute mp3 energy or S+D part of mp4 energy
  std::shared_ptr<PSIO> psio(new PSIO());

  if (iter == 1){
     if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempt[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempt;
     }
     emp3_os = 0.0;
     emp3_ss = 0.0;
     for (a = 0; a < v; a++) {
         for (b = 0; b < v; b++) {
             for (i = 0; i < o; i++) {
                 for (j = 0; j < o; j++) {
                     emp3_os += tb[a*o*o*v+b*o*o+i*o+j] * tempv[a*o*o*v+b*o*o+i*o+j];
                     emp3_ss += (tb[a*o*o*v+b*o*o+i*o+j] - tb[a*o*o*v+b*o*o+j*o+i]) * tempv[a*o*o*v+b*o*o+i*o+j];
                 }
             }
         }
     }
     emp3 = emp3_os + emp3_ss;
  }else if (iter == 2){
     if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempt;
     }
     emp4_sd_os = 0.0;
     emp4_sd_ss = 0.0;
     for (a = 0; a < v; a++) {
         for (b = 0; b < v; b++) {
             for (i = 0; i < o; i++) {
                 for (j = 0; j < o; j++) {
                     emp4_sd_os += tb[a*o*o*v+b*o*o+i*o+j] * tempv[a*o*o*v+b*o*o+i*o+j];
                     emp4_sd_ss += (tb[a*o*o*v+b*o*o+i*o+j] - tb[a*o*o*v+b*o*o+j*o+i]) * tempv[a*o*o*v+b*o*o+i*o+j];
                 }
             }
         }
     }
     emp4_sd = emp4_sd_os + emp4_sd_ss;
  }else if (iter == 3){
     if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempt;
     }
     emp4_q_os = 0.0;
     emp4_q_ss = 0.0;
     for (a = 0; a < v; a++) {
         for (b = 0; b < v; b++) {
             for (i = 0; i < o; i++) {
                 for (j = 0; j < o; j++) {
                     emp4_q_os += tb[a*o*o*v+b*o*o+i*o+j] * tempv[a*o*o*v+b*o*o+i*o+j];
                     emp4_q_ss += (tb[a*o*o*v+b*o*o+i*o+j] - tb[a*o*o*v+b*o*o+j*o+i]) * tempv[a*o*o*v+b*o*o+i*o+j];
                 }
             }
         }
     }
     emp4_q = emp4_q_os + emp4_q_ss;
  }

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (iter == 0) {

      for (i=0; i<o; i++){
          di = - eps[i];
          for (j=0; j<o; j++){
              dij = di-eps[j];
              for (a=o; a<rs; a++){
                  dija = dij + eps[a];
                  for (b=o; b<rs; b++){
                      dijab = dija + eps[b];
                      iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                      ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                      tnew = - integrals[iajb]/dijab;
                      tempt[ijab] = tnew;

                  }
              }
          }
      }
      emp2_os = 0.0;
      emp2_ss = 0.0;
      for (i=0; i<o; i++){
          for (j=0; j<o; j++){
              for (a=o; a<rs; a++){
                  for (b=o; b<rs; b++){
                      iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                      ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                      tnew = - integrals[iajb]/dijab;
                      emp2_os += tempt[ijab] * integrals[iajb];
                      emp2_ss += (tempt[ijab] - tempt[(a-o)*o*o*v+(b-o)*o*o+j*o+i]) * integrals[iajb];

                  }
              }
          }
      }
      emp2 = emp2_os + emp2_ss;

      psio->open(PSIF_DCC_T2,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
      psio->write_entry(PSIF_DCC_T2,"first",(char*)&tempt[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);

      if (!t2_on_disk) C_DCOPY(o*o*v*v,tempt,1,tb,1);

  }
  else if (iter == 1) {
      for (i=0; i<o; i++){
          di = - eps[i];
          for (j=0; j<o; j++){
              dij = di-eps[j];
              for (a=o; a<rs; a++){
                  dija = dij + eps[a];
                  for (b=o; b<rs; b++){
                      dijab = dija + eps[b];
                      iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                      ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                      tnew = - tempv[ijab]/dijab;
                      tempt[ijab] = tnew;
                  }
              }
          }
      }
  }
  psio.reset();
}

/**
 *  Build and use I2ijkl
 */
void CoupledCluster::I2ijkl_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
  }

  psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJKL,"E2ijkl",(char*)&integrals[0],o*o*o*o*sizeof(double));
  psio->close(PSIF_DCC_IJKL,1);

  F_DGEMM('n','n',o*o,v*v,o*o,0.5,integrals,o*o,tempt,o*o,0.0,tempv,o*o);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+b*v*o*o+a*o*o+i,o,tempt+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();

}
/**
 *  Build and use I2'iajk
 */
void CoupledCluster::I2piajk_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK2,"E2ijak2",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_DCC_IJAK2,1);

  F_DGEMM('n','n',o*o*v,v,o,-1.0,tempv,o*o*v,t1,o,0.0,tempt,o*o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  C_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempt+b*v*o*o+a*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();
}
/**
 *  Use Vabcd1
 */
void CoupledCluster::Vabcd1_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=a+1; b<v; b++){
                  tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                     tempt[a*o*o*v+b*o*o+i*o+j]+tempt[b*o*o*v+a*o*o+i*o+j];
              }
              tempv[Position(a,a)*o*(o+1)/2+Position(i,j)] =
                 tempt[a*o*o*v+a*o*o+i*o+j];
          }
      }
  }
  psio->open(PSIF_DCC_ABCD1,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_DCC_ABCD1,"E2abcd1",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j=ntiles-1;
  psio->read(PSIF_DCC_ABCD1,"E2abcd1",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_DCC_ABCD1,1);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();

}
/**
 *  Use Vabcd2
 */
void CoupledCluster::Vabcd2_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  int sg,sg2;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=a; b<v; b++){
                  tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                    tempt[a*o*o*v+b*o*o+i*o+j]-tempt[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }
  psio->open(PSIF_DCC_ABCD2,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_DCC_ABCD2,"E2abcd2",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j = ntiles-1;
  psio->read(PSIF_DCC_ABCD2,"E2abcd2",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_DCC_ABCD2,1);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          if (a>b) sg2 = -1;
          else     sg2 = 1;
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  if (i>j) sg = -1;
                  else     sg = 1;
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*sg2*sg*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  //psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();
}
/**
 *  Build and use I2iabj
 */
void CoupledCluster::I2iabj_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  C_DCOPY(o*o*v*v,integrals,1,tempv,1);

  // use I2iabj
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }
  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  tempt[id++] = 2*tb[a*o*o*v+b*o*o+i*o+j]-tb[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  // if we KNOW this is the first diagram, we don't need to read in the old
  // residual.
  //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] = integrals[j*o*v*v+b*v*o+i*v+a] + integrals[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();

}
/**
 *  Build and use I2iajb
 */
void CoupledCluster::I2iajb_linear(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);

  // use I2iajb
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  integrals[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,integrals,o*v,0.0,tempv,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  integrals[id++] += tempv[j*o*v*v+b*v*o+i*v+a] + tempv[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // use I2iajb
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }

  for (j=0,id=0; j<o; j++){
      for (a=0; a<v; a++){
          for (i=0; i<o; i++){
              for (b=0; b<v; b++){
                  tempv[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,tempv,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (i=0; i<o; i++){
                  tempt[id++] += integrals[j*o*v*v+b*v*o+i*v+a] + integrals[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();
}

/**
 *  Tasks: linear diagrams
 */
void CoupledCluster::DefineLinearTasks(){
  LTasklist = new CCTask[1000];
  LParams   = new CCTaskParams[1000];
  long int o = ndoccact;
  long int v = nvirt;

  nltasks=0;

  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::I2iabj_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::I2iajb_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::I2ijkl_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::I2piajk_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmeni_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmaef_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::CPU_I2p_abci_refactored_term1_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmeai_linear;
  LTasklist[nltasks++].func  = &psi::fnocc::CoupledCluster::Vabcd1_linear;
  // this is the last diagram that contributes to doubles residual,
  // so we can keep it in memory rather than writing and rereading
  LTasklist[nltasks++].func        = &psi::fnocc::CoupledCluster::Vabcd2_linear;
}

/**
 *  Tasks: quadratic diagrams
 */
void CoupledCluster::DefineQuadraticTasks(){
  QTasklist = new CCTask[1000];
  QParams   = new CCTaskParams[1000];

  nqtasks=0;

  QTasklist[nqtasks++].func  = &psi::fnocc::CoupledCluster::I2iabj_quadratic;
  QTasklist[nqtasks++].func  = &psi::fnocc::CoupledCluster::I2ijkl_quadratic;
  QTasklist[nqtasks++].func  = &psi::fnocc::CoupledCluster::CPU_I1ab_quadratic;
  QTasklist[nqtasks++].func  = &psi::fnocc::CoupledCluster::CPU_I1pij_I1ia_lessmem_quadratic;
  QTasklist[nqtasks++].func  = &psi::fnocc::CoupledCluster::I2iajb_quadratic;

}

}} // end of namespace psi
