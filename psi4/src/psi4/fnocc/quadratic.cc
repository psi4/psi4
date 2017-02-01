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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include"psi4/libqt/qt.h"
#include<sys/times.h>
#include "psi4/libpsio/psio.hpp"
#include "psi4/psifiles.h"
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"blas.h"
#include"ccsd.h"

using namespace psi;

namespace psi{ namespace fnocc{

/**
  * Build and use I(a,b)
  */
void CoupledCluster::CPU_I1ab_quadratic(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int b,m,n,e,a,id=0;
  // build I1(a,b)
  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              C_DCOPY(v,tb+e*v*o*o+m*o+n,o*o,tempt+m*o*v*v+e*o*v+n*v,1);
          }
      }
  }
  C_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              C_DAXPY(v,-0.5,integrals+m*o*v*v+n*v+e,o*v,tempv+m*o*v*v+e*o*v+n*v,1);
          }
      }
  }
  F_DGEMM('n','t',v,v,o*o*v,-2.0,tempv,v,tempt,v,0.0,I1,v);

  long int i,j,l,k,c,d;

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (l=0,id=0; l<o; l++){
      for (c=0; c<v; c++){
          for (k=0; k<o; k++){
              C_DCOPY(v,tb+c*o*o+l*o+k,v*o*o,tempt+l*o*v*v+c*o*v+k*v,1);
          }
      }
  }
  // use I1(a,b) for doubles residual:
  F_DGEMM('t','n',v,o*o*v,v,1.0,I1,v,tempt,v,0.0,tempv,v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+a*v*o+i*v+b,v*v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+i*v*v*o+b*v*o+a,v,tempt+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

/**
  * Build and use I(i,j), I'(i,j), and I(i,a)
  */
void CoupledCluster::CPU_I1pij_I1ia_lessmem_quadratic(CCTaskParams params){

  long int o = ndoccact;
  long int v = nvirt;
  long int m,j,e,f,i,a,b;//,one=1;
  long int ov2 = o*v*v;
  long int id=0;
  std::shared_ptr<PSIO> psio(new PSIO());

  // no singles
  // build I1(i,a). n^4
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  //C_DCOPY(o*o*v*v,integrals,1,tempv,1);
  //for (i=0; i<o; i++){
  //    for (a=0; a<v; a++){
  //        for (m=0; m<o; m++){
  //            C_DAXPY(v,-0.5,integrals+i*o*v*v+m*v+a,o*v,tempv+i*v*v*o+a*v*o+m*v,1);
  //        }
  //    }
  //}
  //for (i=0; i<o; i++) C_DCOPY(v,t1+i,o,tempt+i*v,1);
  //F_DGEMV('t',o*v,o*v,2.0,tempv,o*v,tempt,1,0.0,I1,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  memset((void*)tempt,'\0',o*o*v*v);
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
              C_DAXPY(v,-0.5,tb+e*o*o*v+j*o+m,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  // no singles
  // use I1(i,a) -> w1
  //F_DGEMV('n',o*v,o*v,2.0,tempt,o*v,I1,1,0.0,tempv,1);
  //for (i=0; i<o; i++){
  //    C_DAXPY(v,1.0,tempv+i*v,1,w1+i,o);
  //}

  // build I1'(i,j)
  F_DGEMM('t','n',o,o,ov2,2.0,tempt,ov2,integrals,ov2,0.0,I1p,o);

  // no singles
  // only n^4
  //psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&tempt[0],o*o*o*v*sizeof(double));
  //psio->close(PSIF_IJAK,1);
  //id=0;
  //for (i=0; i<o; i++){
  //    for (j=0; j<o; j++){
  //        for (e=0; e<v; e++){
  //            C_DCOPY(o,tempt+i*o*v+j*v+e,o*o*v,tempv+i*o*o*v+j*o*v+e*o,1);
  //            C_DAXPY(o,-2.0,tempt+i*o*o*v+j*v+e,o*v,tempv+i*o*o*v+j*o*v+e*o,1);
  //        }
  //    }
  //}
  //F_DGEMV('t',o*v,o*o,-1.0,tempv,o*v,t1,1,1.0,I1p,1);

  // no singles
  // use I1'(i,j) for singles residual. (n^3)
  //F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);
  //F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);

  // no singles
  // build I1(i,j)
  //F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);
  //F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  F_DGEMM('n','t',o,ov2,o,-1.0,I1p,o,tempt,ov2,0.0,tempv,o);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+a*o*o*v+b*o+i,v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+b*o*o*v+i*v*o+a*o,1,tempt+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

/**
 *  Build and use I2ijkl
 */
void CoupledCluster::I2ijkl_quadratic(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
  }

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  for (j=0; j<o; j++){
      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              C_DCOPY(v,integrals+j*o*v*v+b*o*v+i*v,1,tempv+j*o*v*v+i*v*v+b*v,1);
          }
      }
  }

  // no linear terms
  //psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_IJKL,"E2ijkl",(char*)&integrals[0],o*o*o*o*sizeof(double));
  //psio->close(PSIF_DCC_IJKL,1);
  //F_DGEMM('n','n',o*o,o*o,v*v,1.0,tempt,o*o,tempv,v*v,1.0,integrals,o*o);
  F_DGEMM('n','n',o*o,o*o,v*v,1.0,tempt,o*o,tempv,v*v,0.0,integrals,o*o);

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
 *  Build and use I2iabj
 */
void CoupledCluster::I2iabj_quadratic(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tb+b*v*o*o+j*o+i,o*o,tempt+i*o*v*v+b*o*v+j*v,1);
          }
      }
  }

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  // no linear terms!
  //C_DCOPY(o*o*v*v,integrals,1,tempv,1);
  //F_DGEMM('n','n',o*v,o*v,o*v,-0.5,tempt,o*v,integrals,o*v,1.0,tempv,o*v);
  F_DGEMM('n','n',o*v,o*v,o*v,-0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);

  // contribute to intermediate
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_DCC_TEMP,"temporary",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  C_DCOPY(o*o*v*v,tempt,1,tempv,1);
  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              C_DAXPY(v,-0.5,tempt+i*v*v*o+j*v+b,v*o,tempv+i*o*v*v+b*o*v+j*v,1);
          }
      }
  }

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempt;
  }

  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tb+a*o*o+j*o+i,v*o*o,integrals+i*v*v*o+a*v*o+j*v,1);
          }
      }
  }
  F_DGEMM('n','n',o*v,o*v,o*v,1.0,integrals,o*v,tempv,o*v,0.0,tempt,o*v);

  // contribute to intermediate
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_TEMP,"temporary",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,0);
  C_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);

  // use I2iabj
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }
  for (j=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DCOPY(v,tb+b*o*o+i*o+j,o*o*v,tempt+j*o*v*v+b*o*v+i*v,1);
              C_DAXPY(v,-0.5,tb+b*o*o*v+i*o+j,o*o,tempt+j*o*v*v+b*o*v+i*v,1);
          }
      }
  }
  F_DGEMM('n','n',o*v,o*v,o*v,2.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);

  // contribute to residual
  // this is the first diagram, we don't need to read in the old residual.
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));

  memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,integrals+b*v*o+i*v+a,o*v*v,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,integrals+i*o*v*v+a*v*o+b,v,tempt+a*o*o*v+b*o*o+i*o,1);
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
void CoupledCluster::I2iajb_quadratic(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tb+b*o*o*v+j*o+i,o*o,integrals+i*o*v*v+b*o*v+j*v,1);
          }
      }
  }
  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              C_DCOPY(v,tempt+i*v*v*o+j*v+b,o*v,tempv+i*o*v*v+b*o*v+j*v,1);
          }
      }
  }

  // no linear terms!
  //psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempt[0],o*o*v*v*sizeof(double));
  //psio->close(PSIF_DCC_IJAB,1);
  //F_DGEMM('n','n',o*v,o*v,o*v,-0.5,integrals,o*v,tempv,o*v,1.0,tempt,o*v);
  F_DGEMM('n','n',o*v,o*v,o*v,-0.5,integrals,o*v,tempv,o*v,0.0,tempt,o*v);

  // use I2iajb
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (j=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DCOPY(v,tb+b*v*o*o+j*o+i,o*o,integrals+j*o*v*v+b*o*v+i*v,1);
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,integrals,o*v,0.0,tempv,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+b*v*o+i*v+a,o*v*v,integrals+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+i*o*v*v+a*v*o+b,v,integrals+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // use I2iajb

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }

  for (j=0; j<o; j++){
      for (a=0; a<v; a++){
          for (i=0; i<o; i++){
              C_DCOPY(v,tb+a*o*o+j*o+i,o*o*v,tempv+j*o*v*v+a*o*v+i*v,1);
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,tempv,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));

  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              C_DAXPY(o,1.0,integrals+j*o*v*v+b*v*o+a,v,tempv+a*o*o*v+b*o*o+j*o,1);
              C_DAXPY(o,1.0,integrals+a*v*o+j*v+b,o*v*v,tempv+a*o*o*v+b*o*o+j*o,1);
          }
      }
  }

  // last diagram.  put in tempv - no need to write
  //psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

}} // end of namespaces
