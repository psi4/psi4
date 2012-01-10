#include"psi4-dec.h"
#include"ccsd.h"
#include <libplugin/plugin.h>
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include"blas.h"
#ifdef _OPENMP
#include<omp.h>
#endif

using namespace psi;
using namespace boost;

namespace psi{
PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options);

PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options){

  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                      CCSD(T)                        *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n");
  fflush(outfile);

  int o = ccsd->ndoccact;
  int v = ccsd->nvirt_no;

  double *t1 = ccsd->t1;
  double *F  = ccsd->eps;
  double *E2ijak,**E2abci;
  E2ijak = (double*)malloc(o*o*o*v*sizeof(double));
  int nthreads = 1;
  #ifdef _OPENMP
      nthreads = omp_get_max_threads();
  #endif

  if (options["NUM_THREADS"].has_changed())
     nthreads = options.get_int("NUM_THREADS");

  long int memory = Process::environment.get_memory();
  memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+3L*nthreads*v*v*v);

  fprintf(outfile,"        num_threads =             %9i\n",nthreads);
  fprintf(outfile,"        available memory =     %9.2lf mb\n",Process::environment.get_memory()/1024./1024.);
  fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
           8.*(2.*o*o*v*v+1.*o*o*o*v+(3.*nthreads)*v*v*v+1.*o*v)/1024./1024.);
  fprintf(outfile,"\n");
  fflush(outfile);

  if (memory<0){
     while (memory<0 && nthreads>0){
           memory += 8L*3L*v*v*v;
           nthreads--;
     }
     if (nthreads<1){
        fprintf(outfile,"        Error: not enough memory.\n");
        fprintf(outfile,"\n");
        fprintf(outfile,"        Setting num_threads = 1 will reduce required memory to %7.2lf mb\n",
             8.*(2.*o*o*v*v+1.*o*o*o*v+3.*v*v*v+1.*o*v)/1024./1024.);
        fprintf(outfile,"\n");
        fflush(outfile);
        return Failure;
     }
     fprintf(outfile,"        Not enough memory.  Decreasing num_threads ... \n");
     fprintf(outfile,"\n");
     fprintf(outfile,"        num_threads =             %9i\n",nthreads);
     fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
              8.*(2.*o*o*v*v+1.*o*o*o*v+(3.*nthreads)*v*v*v+1.*o*v)/1024./1024.);
     fprintf(outfile,"\n");
     fflush(outfile);
  }

  int nijk = 0;
  for (int i=0; i<o; i++){
      for (int j=0; j<=i; j++){
          for (int k=0; k<=j; k++){
              nijk++;
          }
      }
  }
  int**ijk = (int**)malloc(nijk*sizeof(int*));
  nijk = 0;
  for (int i=0; i<o; i++){
      for (int j=0; j<=i; j++){
          for (int k=0; k<=j; k++){
              ijk[nijk] = (int*)malloc(3*sizeof(int));
              ijk[nijk][0] = i;
              ijk[nijk][1] = j;
              ijk[nijk][2] = k;
              nijk++;
          }
      }
  }
  fprintf(outfile,"        Number of ijk combinations: %i\n",nijk);
  fprintf(outfile,"\n");
  fflush(outfile);
  

  E2abci = (double**)malloc(nthreads*sizeof(double*));
  // some v^3 intermediates
  double **Z  = (double**)malloc(nthreads*sizeof(double*));
  double **Z2 = (double**)malloc(nthreads*sizeof(double*));
  //double **Y  = (double**)malloc(nthreads*sizeof(double*));
  //double **Z3 = (double**)malloc(nthreads*sizeof(double*));

  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(v*v*v*sizeof(double));
      Z[i]      = (double*)malloc(v*v*v*sizeof(double));
      Z2[i]     = (double*)malloc(v*v*v*sizeof(double));
      //Y[i]      = (double*)malloc(v*v*v*sizeof(double));
      //Z3[i]     = (double*)malloc(v*v*v*sizeof(double));
  }

  boost::shared_ptr<PSIO> psio(new PSIO());

  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&E2ijak[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);

  double *tempt = (double*)malloc(o*o*v*v*sizeof(double));

  if (ccsd->t2_on_disk){
     ccsd->tb = (double*)malloc(o*o*v*v*sizeof(double));
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&ccsd->tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }

  for (int a=0; a<v*v; a++){
      F_DCOPY(o*o,ccsd->tb+a*o*o,1,tempt+a,v*v);
  }

  // might as well use t2's memory
  double*E2klcd = ccsd->tb;
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd", (char*)&E2klcd[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  double *renorm = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;
  for (int i=0; i<nthreads; i++) renorm[i] = 0.0;
  fprintf(outfile,"        Computing (T) correction... \n");
  fflush(outfile);
  psio->open(PSIF_ABCI,PSIO_OPEN_OLD);

  

  //ccsd->helper_->num_gpus = 0;
  #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
  for (int ind=0; ind<nijk; ind++){
      int i = ijk[ind][0];
      int j = ijk[ind][1];
      int k = ijk[ind][2];

      int thread = 0;
      #ifdef _OPENMP
          thread = omp_get_thread_num();
      #endif

      psio_address addr;
      #pragma omp critical
      {
          addr = psio_get_address(PSIO_ZERO,(long int)k*v*v*v*sizeof(double));
          psio->read(PSIF_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      }
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+j*v*v*o+i*v*v,v,0.0,Z[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+j*o*o*v+k*o*v,v,tempt+i*v*v*o,v*v,1.0,Z[thread],v,thread);

      //(ab)(ij)
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+i*v*v*o+j*v*v,v,0.0,Z2[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+i*o*o*v+k*o*v,v,tempt+j*v*v*o,v*v,1.0,Z2[thread],v,thread);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v*v+a*v,1,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(bc)(jk)
      #pragma omp critical
      {
          addr = psio_get_address(PSIO_ZERO,(long int)j*v*v*v*sizeof(double));
          psio->read(PSIF_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      }
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+k*v*v*o+i*v*v,v,0.0,Z2[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+k*o*o*v+j*o*v,v,tempt+i*v*v*o,v*v,1.0,Z2[thread],v,thread);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+a*v*v+b,v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ikj)(acb)
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+i*v*v*o+k*v*v,v,0.0,Z2[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+i*o*o*v+j*o*v,v,tempt+k*v*v*o,v*v,1.0,Z2[thread],v,thread);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+a*v+b,v*v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ac)(ik)
      #pragma omp critical
      {
          addr = psio_get_address(PSIO_ZERO,(long int)i*v*v*v*sizeof(double));
          psio->read(PSIF_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      }
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+j*v*v*o+k*v*v,v,0.0,Z2[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+j*o*o*v+i*o*v,v,tempt+k*v*v*o,v*v,1.0,Z2[thread],v,thread);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v+a,v*v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ijk)(abc)
      ccsd->helper_->GPUTiledDGEMM_NoThread('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+k*v*v*o+j*v*v,v,0.0,Z2[thread],v*v,thread);
      ccsd->helper_->GPUTiledDGEMM_NoThread('n','t',v,v*v,o,-1.0,E2ijak+k*o*o*v+i*o*v,v,tempt+j*v*v*o,v*v,1.0,Z2[thread],v,thread);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v*v+a,v,Z[thread]+a*v*v+b*v,1);
          }
      }

      F_DCOPY(v*v*v,Z[thread],1,Z2[thread],1);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  Z2[thread][a*v*v+b*v+c] += (t1[a*o+i]*E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                              t1[b*o+j]*E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                              t1[c*o+k]*E2klcd[i*v*v*o+a*v*o+j*v+b]);
                  Z2[thread][a*v*v+b*v+c] /= (1 + (a==b) + (b==c) + (a==c));
              }
          }
      }

      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bac = b*v*v+a*v+c;
                  long int acb = a*v*v+c*v+b;
                  long int bca = b*v*v+c*v+a;
                  long int cab = c*v*v+a*v+b;
                  long int cba = c*v*v+b*v+a;

                  //Y[thread][abc]  = Z2[thread][abc] + Z2[thread][bca] + Z2[thread][cab];
                  E2abci[thread][abc] = Z2[thread][acb] + Z2[thread][bac] + Z2[thread][cba];
              }
          }
      }
      double dijk = F[i]+F[j]+F[k];
      // separate out these bits to save v^3 storage
      for (int a=0; a<v; a++){
          double dijka = dijk-F[a+o];
          for (int b=0; b<=a; b++){
              double dijkab = dijka-F[b+o];
              for (int c=0; c<=b; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bca = b*v*v+c*v+a;
                  long int cab = c*v*v+a*v+b;
                  long int acb = a*v*v+c*v+b;
                  long int bac = b*v*v+a*v+c;
                  long int cba = c*v*v+b*v+a;
                  double dum      = Z[thread][abc]*Z2[thread][abc] + Z[thread][acb]*Z2[thread][acb]
                                  + Z[thread][bac]*Z2[thread][bac] + Z[thread][bca]*Z2[thread][bca]
                                  + Z[thread][cab]*Z2[thread][cab] + Z[thread][cba]*Z2[thread][cba];

                  dum            =  (E2abci[thread][abc])
                                 * ((Z[thread][abc] + Z[thread][bca] + Z[thread][cab])*-2.0
                                 +  (Z[thread][acb] + Z[thread][bac] + Z[thread][cba]))
                                 + 3.0*dum;
                  //dum            = (Y[thread][abc] - 2.0*E2abci[thread][abc])
                  //               * (Z[thread][abc] + Z[thread][bca] + Z[thread][cab])
                  //               + (E2abci[thread][abc] - 2.0*Y[thread][abc])
                  //               * (Z[thread][acb] + Z[thread][bac] + Z[thread][cba])
                  //               + 3.0*dum;
                  double denom = dijkab-F[c+o];
                  etrip[thread] += dum/denom*( 2-((i==j)+(j==k)+(i==k)) );
              }
          }
      }
      // the second bit
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bac = b*v*v+a*v+c;
                  long int acb = a*v*v+c*v+b;
                  long int bca = b*v*v+c*v+a;
                  long int cab = c*v*v+a*v+b;
                  long int cba = c*v*v+b*v+a;

                  E2abci[thread][abc]  = Z2[thread][abc] + Z2[thread][bca] + Z2[thread][cab];
              }
          }
      }
      for (int a=0; a<v; a++){
          double dijka = dijk-F[a+o];
          for (int b=0; b<=a; b++){
              double dijkab = dijka-F[b+o];
              for (int c=0; c<=b; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bca = b*v*v+c*v+a;
                  long int cab = c*v*v+a*v+b;
                  long int acb = a*v*v+c*v+b;
                  long int bac = b*v*v+a*v+c;
                  long int cba = c*v*v+b*v+a;

                  double dum     = (E2abci[thread][abc])
                                 * (Z[thread][abc] + Z[thread][bca] + Z[thread][cab]
                                 + (Z[thread][acb] + Z[thread][bac] + Z[thread][cba])*-2.0);

                  double denom = dijkab-F[c+o];
                  etrip[thread] += dum/denom*( 2-((i==j)+(j==k)+(i==k)) );
              }
          }
      }
  }


  double et = 0.0;
  for (int i=0; i<nthreads; i++) et += etrip[i];

  psio->close(PSIF_ABCI,1);
  fprintf(outfile,"\n");
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"        (T) energy                   %20.12lf\n",et);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        (T) energy                   %20.12lf %20.12lf\n",et,et*ccsd->scale_t);
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"        MP2 correlation energy       %20.12lf\n",ccsd->emp2);
  fprintf(outfile,"        CCSD correlation energy      %20.12lf\n",ccsd->eccsd);
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf\n",ccsd->eccsd+et);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf %20.12lf\n",ccsd->eccsd+et,ccsd->eccsd+et*ccsd->scale_t);
  }
  fflush(outfile);
  ccsd->et = et;

  // free memory:
  if (ccsd->t2_on_disk){
     free(ccsd->tb);
  }
  free(E2ijak);
  for (int i=0; i<nthreads; i++){  
      free(E2abci[i]);
      //free(Y[i]);
      free(Z[i]);
      free(Z2[i]);
      //free(Z3[i]);
  }
  //free(Y);
  free(Z);
  free(Z2);
  //free(Z3);
  free(E2abci);
  free(etrip);
  free(renorm);
            
  return Success;
}


} // end of namespace



