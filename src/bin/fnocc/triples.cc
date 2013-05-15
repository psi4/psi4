/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include"ccsd.h"
#include"blas.h"
#include<libmints/wavefunction.h>
#ifdef _OPENMP
   #include<omp.h>
#endif

using namespace psi;

namespace psi{namespace fnocc{

PsiReturnType CoupledCluster::triples(){

  char*name = new char[10];
  char*space = new char[10];
  double fac;
  if (ccmethod == 0) {
     sprintf(name,"CCSD");
     sprintf(space," ");
     fac = 1.0;
  }else if (ccmethod == 1) {
     sprintf(name,"QCISD");
     sprintf(space,"  ");
     fac = 2.0;
  }else{
     sprintf(name,"MP4");
     sprintf(space,"");
     fac = 0.0;
  }

  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                  %8s(T)                        *\n",name);
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n");
  fflush(outfile);

  int o = ndoccact;
  int v = nvirt_no;

  double *F  = eps;
  double *E2ijak,**E2abci;
  E2ijak = (double*)malloc(o*o*o*v*sizeof(double));
  int nthreads = 1;
  #ifdef _OPENMP
      nthreads = omp_get_max_threads();
  #endif

  long int memory = Process::environment.get_memory();
  if (options_["MEMORY"].has_changed()){
     memory  = options_.get_int("MEMORY");
     memory *= (long int)1024*1024;
  }
  memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+3L*nthreads*v*v*v);

  fprintf(outfile,"        num_threads =             %9i\n",nthreads);
  fprintf(outfile,"        available memory =     %9.2lf mb\n",memory/1024./1024.);
  fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
           8.*(2.*o*o*v*v+1.*o*o*o*v+(3.*nthreads)*v*v*v+1.*o*v)/1024./1024.);
  fprintf(outfile,"\n");
  fflush(outfile);

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

  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(v*v*v*sizeof(double));
      Z[i]      = (double*)malloc(v*v*v*sizeof(double));
      Z2[i]     = (double*)malloc(v*v*v*sizeof(double));
  }

  boost::shared_ptr<PSIO> psio(new PSIO());

  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&E2ijak[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_DCC_IJAK,1);

  double *tempt = (double*)malloc(o*o*v*v*sizeof(double));

  // first-order amplitudes for mp4
  if (ccmethod == 2) {
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }

  for (int a=0; a<v*v; a++){
      F_DCOPY(o*o,tb+a*o*o,1,tempt+a,v*v);
  }

  // might as well use t2's memory
  double*E2klcd = tb;
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb", (char*)&E2klcd[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;
  fprintf(outfile,"        Computing (T) correction...\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"        %% complete  total time\n");
  fflush(outfile);

  time_t stop,start = time(NULL);
  int pct10,pct20,pct30,pct40,pct50,pct60,pct70,pct80,pct90;
  pct10=pct20=pct30=pct40=pct50=pct60=pct70=pct80=pct90=0;

  /**
    *  if there is enough memory to explicitly thread, do so
    */
  #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
  for (int ind=0; ind<nijk; ind++){
      int i = ijk[ind][0];
      int j = ijk[ind][1];
      int k = ijk[ind][2];

      int thread = 0;
      #ifdef _OPENMP
          thread = omp_get_thread_num();
      #endif

      boost::shared_ptr<PSIO> mypsio(new PSIO());
      mypsio->open(PSIF_DCC_ABCI,PSIO_OPEN_OLD);

      psio_address addr = psio_get_address(PSIO_ZERO,(long int)k*v*v*v*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+j*v*v*o+i*v*v,v,0.0,Z[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+j*o*o*v+k*o*v,v,tempt+i*v*v*o,v*v,1.0,Z[thread],v);

      //(ab)(ij)
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+i*v*v*o+j*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+i*o*o*v+k*o*v,v,tempt+j*v*v*o,v*v,1.0,Z2[thread],v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v*v+a*v,1,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(bc)(jk)
      addr = psio_get_address(PSIO_ZERO,(long int)j*v*v*v*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+k*v*v*o+i*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+k*o*o*v+j*o*v,v,tempt+i*v*v*o,v*v,1.0,Z2[thread],v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+a*v*v+b,v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ikj)(acb)
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+i*v*v*o+k*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+i*o*o*v+j*o*v,v,tempt+k*v*v*o,v*v,1.0,Z2[thread],v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+a*v+b,v*v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ac)(ik)
      addr = psio_get_address(PSIO_ZERO,(long int)i*v*v*v*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],v*v*v*sizeof(double),addr,&addr);
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+j*v*v*o+k*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+j*o*o*v+i*o*v,v,tempt+k*v*v*o,v*v,1.0,Z2[thread],v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v+a,v*v,Z[thread]+a*v*v+b*v,1);
          }
      }

      //(ijk)(abc)
      F_DGEMM('t','t',v*v,v,v,1.0,E2abci[thread],v,tempt+k*v*v*o+j*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,v*v,o,-1.0,E2ijak+k*o*o*v+i*o*v,v,tempt+j*v*v*o,v*v,1.0,Z2[thread],v);
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              F_DAXPY(v,1.0,Z2[thread]+b*v*v+a,v,Z[thread]+a*v*v+b*v,1);
          }
      }

      F_DCOPY(v*v*v,Z[thread],1,Z2[thread],1);
      for (int a=0; a<v; a++){
          double tai = t1[a*o+i];
          for (int b=0; b<v; b++){
              int ab = 1+(a==b);
              double tbj = t1[b*o+j];
              double E2iajb = E2klcd[i*v*v*o+a*v*o+j*v+b];
              for (int c=0; c<v; c++){
                  Z2[thread][a*v*v+b*v+c] += fac*(tai      *E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                              tbj      *E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                              t1[c*o+k]*E2iajb);
                  Z2[thread][a*v*v+b*v+c] /= (ab + (b==c) + (a==c));
              }
          }
      }

      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bac = b*v*v+a*v+c;
                  long int acb = a*v*v+c*v+b;
                  long int cba = c*v*v+b*v+a;

                  E2abci[thread][abc] = Z2[thread][acb] + Z2[thread][bac] + Z2[thread][cba];
              }
          }
      }
      double dijk = F[i]+F[j]+F[k];
      int ijkfac = ( 2-((i==j)+(j==k)+(i==k)) );
      // separate out these bits to save v^3 storage
      double tripval = 0.0;
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
                  double denom = dijkab-F[c+o];
                  tripval += dum/denom;
              }
          }
      }
      etrip[thread] += tripval*ijkfac;
      // the second bit
      for (int a=0; a<v; a++){
          for (int b=0; b<v; b++){
              for (int c=0; c<v; c++){
                  long int abc = a*v*v+b*v+c;
                  long int bca = b*v*v+c*v+a;
                  long int cab = c*v*v+a*v+b;

                  E2abci[thread][abc]  = Z2[thread][abc] + Z2[thread][bca] + Z2[thread][cab];
              }
          }
      }
      tripval = 0.0;
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
                  tripval += dum/denom;
              }
          }
      }
      etrip[thread] += tripval*ijkfac;
      // print out update 
      if (thread==0){
         int print = 0;
         stop = time(NULL);
         if ((double)ind/nijk >= 0.1 && !pct10){      pct10 = 1; print=1;}
         else if ((double)ind/nijk >= 0.2 && !pct20){ pct20 = 1; print=1;}
         else if ((double)ind/nijk >= 0.3 && !pct30){ pct30 = 1; print=1;}
         else if ((double)ind/nijk >= 0.4 && !pct40){ pct40 = 1; print=1;}
         else if ((double)ind/nijk >= 0.5 && !pct50){ pct50 = 1; print=1;}
         else if ((double)ind/nijk >= 0.6 && !pct60){ pct60 = 1; print=1;}
         else if ((double)ind/nijk >= 0.7 && !pct70){ pct70 = 1; print=1;}
         else if ((double)ind/nijk >= 0.8 && !pct80){ pct80 = 1; print=1;}
         else if ((double)ind/nijk >= 0.9 && !pct90){ pct90 = 1; print=1;}
         if (print){
            fprintf(outfile,"              %3.1lf  %8d s\n",100.0*ind/nijk,(int)stop-(int)start);
            fflush(outfile);
         }
      }
      mypsio->close(PSIF_DCC_ABCI,1);
      mypsio.reset();
  }

  double myet = 0.0;
  for (int i=0; i<nthreads; i++) myet += etrip[i];

  // ccsd(t) or qcisd(t)
  if (ccmethod <= 1) {
      et = myet;
      fprintf(outfile,"\n");
      fprintf(outfile,"        (T) energy   %s                   %20.12lf\n",space,et);
      fprintf(outfile,"\n");
      fprintf(outfile,"        %s(T) correlation energy       %20.12lf\n",name,eccsd+et);
      fprintf(outfile,"      * %s(T) total energy             %20.12lf\n",name,eccsd+et+escf);
      fprintf(outfile,"\n");
  }else {
      emp4_t = myet;
      fprintf(outfile,"\n");
      fprintf(outfile,"        MP4(T) correlation energy:         %20.12lf\n",emp4_t);
      fprintf(outfile,"\n");
      fprintf(outfile,"        MP4(SDTQ) correlation energy:      %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+emp4_t);
      fprintf(outfile,"      * MP4(SDTQ) total energy:            %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+emp4_t+escf);
      fprintf(outfile,"\n");
  }
  fflush(outfile);

  // free memory:
  free(E2ijak);
  free(tempt);
  for (int i=0; i<nthreads; i++){  
      free(E2abci[i]);
      free(Z[i]);
      free(Z2[i]);
  }
  free(Z);
  free(Z2);
  free(E2abci);
  free(etrip);
  delete name;
  delete space;
            
  return Success;
}

}} // end of namespaces



