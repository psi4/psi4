#include"psi4-dec.h"
#include<libplugin/plugin.h>
#include<boost/shared_ptr.hpp>
#include<lib3index/dftensor.h>
#include<liboptions/liboptions.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
#include<libchkpt/chkpt.h>
#include<libiwl/iwl.h>
#include<libpsio/psio.hpp>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include"ccsd.h"
#include"blas.h"

using namespace psi;
using namespace boost;

namespace psi{

void DensityFittedIntegrals(){
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  int nocc = ref->doccpi()[0];
  int nvir = ref->nmopi()[0]-ref->doccpi()[0];
  int aocc = nocc-ref->frzcpi()[0];
  int avir = nvir-ref->frzvpi()[0];

  boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, ref->molecule(), "RI_BASIS_CC");

  boost::shared_ptr<DFTensor> DF (new DFTensor(ref->basisset(),auxiliary,ref->Ca(),nocc,nvir,aocc,avir,Process::environment.options));
  SharedMatrix Qoo = DF->Qoo();
  SharedMatrix Qov = DF->Qov();
  SharedMatrix Qvv = DF->Qvv();

  double**Qoo_pointer = Qoo->pointer();
  double**Qov_pointer = Qov->pointer();
  double**Qvv_pointer = Qvv->pointer();

  int nQ = auxiliary->nbf();
  int o = aocc;
  int v = avir;
  int dim = o*o*v*v;
  if (v*v*v>dim) dim = v*v*v;
  double* temp1 = (double*)malloc(dim*sizeof(double));
  double* temp2 = (double*)malloc(dim*sizeof(double));
  double* tempq = (double*)malloc(nQ*v*v*sizeof(double));

  psio_address addr,addr2,addr3;
  boost::shared_ptr<PSIO> psio(new PSIO());

  fprintf(outfile,"\n");
  fprintf(outfile,"  Generate density-fitted integrals:\n");
  fprintf(outfile,"\n");

  // (oo|oo)
  fprintf(outfile,"     (oo|oo) block.......");fflush(outfile);
  F_DGEMM('n','t',o*o,o*o,nQ,1.0,&Qoo_pointer[0][0],o*o,&Qoo_pointer[0][0],o*o,0.0,temp1,o*o);
  for (int i=0; i<o; i++){
      for (int j=0; j<o; j++){
          for (int k=0; k<o; k++){
              for (int l=0; l<o; l++){
                  temp2[i*o*o*o+k*o*o+j*o+l] = temp1[i*o*o*o+j*o*o+k*o+l];
              }
          }
      }
  }
  psio->open(PSIF_IJKL,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_IJKL,"E2ijkl",(char*)&temp2[0],o*o*o*o*sizeof(double));
  psio->close(PSIF_IJKL,1);
  fprintf(outfile,"done.\n");fflush(outfile);


  // (oo|ov) 1
  fprintf(outfile,"     (oo|ov) block.......");fflush(outfile);
  F_DGEMM('n','t',o*o,o*v,nQ,1.0,&Qoo_pointer[0][0],o*o,&Qov_pointer[0][0],o*v,0.0,temp1,o*o);
  for (int i=0; i<o; i++){
      for (int j=0; j<o; j++){
          for (int k=0; k<o; k++){
              for (int a=0; a<v; a++){
                  temp2[j*o*o*v+i*o*v+k*v+a] = temp1[i*o*o*v+a*o*o+j*o+k];
              }
          }
      }
  }
  psio->open(PSIF_IJAK,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_IJAK,"E2ijak",(char*)&temp2[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);

  // (oo|ov) 2
  for (int i=0; i<o; i++){
      for (int j=0; j<o; j++){
          for (int k=0; k<o; k++){
              for (int a=0; a<v; a++){
                  temp2[j*o*o*v+a*o*o+k*o+i] = temp1[i*o*o*v+a*o*o+j*o+k];
              }
          }
      }
  }
  psio->open(PSIF_IJAK2,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_IJAK2,"E2ijak2",(char*)&temp2[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK2,1);
  fprintf(outfile,"done.\n");fflush(outfile);

  // (ov|ov)
  fprintf(outfile,"     (ov|ov) block.......");fflush(outfile);
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,&Qov_pointer[0][0],o*v,&Qov_pointer[0][0],o*v,0.0,temp1,o*v);
  psio->open(PSIF_KLCD,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_KLCD,"E2klcd",(char*)&temp1[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  fprintf(outfile,"done.\n");fflush(outfile);

  // (oo|vv)
  fprintf(outfile,"     (oo|vv) block.......");fflush(outfile);
  F_DGEMM('n','t',v*v,o*o,nQ,1.0,&Qvv_pointer[0][0],v*v,&Qoo_pointer[0][0],o*o,0.0,temp1,v*v);
  for (int k=0; k<o; k++){
      for (int c=0; c<v; c++){
          for (int j=0; j<o; j++){
              for (int a=0; a<v; a++){
                  temp2[k*o*v*v+c*o*v+j*v+a] = temp1[k*o*v*v+j*v*v+a*v+c];
              }
          }
      }
  }
  psio->open(PSIF_AKJC2,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_AKJC2,"E2akjc2",(char*)&temp2[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_AKJC2,1);
  fprintf(outfile,"done.\n");fflush(outfile);

  // (ov|vv) 2, 3, and 5
  fprintf(outfile,"     (ov|vv) block.......");fflush(outfile);
  psio->open(PSIF_ABCI2,PSIO_OPEN_NEW);
  psio->open(PSIF_ABCI3,PSIO_OPEN_NEW);
  psio->open(PSIF_ABCI5,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  addr2 = PSIO_ZERO;
  addr3 = PSIO_ZERO;
  for (int a=0; a<v; a++){
      for (int q=0; q<nQ; q++){
          for (int b=0; b<v; b++){
              tempq[q*v+b] = Qvv_pointer[q][a*v+b];
          }
      }
      F_DGEMM('n','t',o*v,v,nQ,1.0,&Qov_pointer[0][0],o*v,tempq,v,0.0,temp1,o*v);
      psio->write(PSIF_ABCI3,"E2abci3",(char*)&temp1[0],o*v*v*sizeof(double),addr,&addr);
      for (int b=0; b<v; b++){
          for (int i=0; i<o; i++){
              for (int c=0; c<v; c++){
                  temp2[b*v*o+i*v+c] = temp1[c*v*o+i*v+b];
              }
          }
      }
      psio->write(PSIF_ABCI5,"E2abci5",(char*)&temp2[0],o*v*v*sizeof(double),addr2,&addr2);

      F_DAXPY(o*v*v,-2.0,temp1,1,temp2,1);
      psio->write(PSIF_ABCI2,"E2abci2",(char*)&temp2[0],o*v*v*sizeof(double),addr3,&addr3);
  }
  psio->close(PSIF_ABCI2,1);
  psio->close(PSIF_ABCI3,1);
  psio->close(PSIF_ABCI5,1);

  // (ov|vv) 1 and 4
  psio->open(PSIF_ABCI,PSIO_OPEN_NEW);
  psio->open(PSIF_ABCI4,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  addr2 = PSIO_ZERO;
  for (int i=0; i<o; i++){
      for (int q=0; q<nQ; q++){
          for (int b=0; b<v; b++){
              tempq[q*v+b] = Qov_pointer[q][i*v+b];
          }
      }
      F_DGEMM('n','t',v,v*v,nQ,1.0,tempq,v,&Qvv_pointer[0][0],v*v,0.0,temp1,v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  temp2[a*v*v+b*v+c] = temp1[a*v*v+c*v+b];
              }
          }
      }
      psio->write(PSIF_ABCI,"E2abci",(char*)&temp2[0],v*v*v*sizeof(double),addr,&addr);
      psio->write(PSIF_ABCI4,"E2abci4",(char*)&temp1[0],v*v*v*sizeof(double),addr2,&addr2);
  }
  psio->close(PSIF_ABCI,1);
  psio->close(PSIF_ABCI4,1);
  fprintf(outfile,"done.\n");fflush(outfile);

  //(vv|vv)
  fprintf(outfile,"     (vv|vv) block.......");fflush(outfile);
  int count = 0;
  addr = PSIO_ZERO;
  addr2 = PSIO_ZERO;
  psio->open(PSIF_ABCD1,PSIO_OPEN_NEW);
  psio->open(PSIF_ABCD2,PSIO_OPEN_NEW);
  for (int a=0; a<v*v; a++){
      for (int q=0; q<nQ; q++){
          tempq[a*nQ+q] = Qvv_pointer[q][a];
      }
  }
  for (int b=0; b<v; b++){
      for (int a=0; a<=b; a++){
          for (int d=0; d<v; d++){
              for (int c=0; c<=d; c++){
                  //double acbd = F_DDOT(nQ,&Qvv_pointer[0][a*v+c],v*v,&Qvv_pointer[0][b*v+d],v*v);
                  //double adbc = F_DDOT(nQ,&Qvv_pointer[0][a*v+d],v*v,&Qvv_pointer[0][b*v+c],v*v);
                  double acbd = F_DDOT(nQ,tempq+a*v*nQ+c*nQ,1,tempq+b*v*nQ+d*nQ,1);
                  double adbc = F_DDOT(nQ,tempq+a*v*nQ+d*nQ,1,tempq+b*v*nQ+c*nQ,1);
                  temp1[count]   = acbd+adbc;
                  temp2[count++] = acbd-adbc;
                  if (count == dim){
                     psio->write(PSIF_ABCD1,"E2abcd1",(char*)&temp1[0],count*sizeof(double),addr,&addr);
                     psio->write(PSIF_ABCD2,"E2abcd2",(char*)&temp2[0],count*sizeof(double),addr2,&addr2);
                     count = 0;
                  }
              }
          }
      }
  }
  if (count!=0){
     psio->write(PSIF_ABCD1,"E2abcd1",(char*)&temp1[0],count*sizeof(double),addr,&addr);
     psio->write(PSIF_ABCD2,"E2abcd2",(char*)&temp2[0],count*sizeof(double),addr2,&addr2);
     count = 0;
  }
  psio->close(PSIF_ABCD1,1);
  psio->close(PSIF_ABCD2,1);
  fprintf(outfile,"done.\n");fflush(outfile);
  fprintf(outfile,"\n");fflush(outfile);

  free(tempq);
  free(temp1);
  free(temp2);
}


} // end of namespace psi
