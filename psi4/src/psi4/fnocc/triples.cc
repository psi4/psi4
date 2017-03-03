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
#include "psi4/libmints/wavefunction.h"
#include"psi4/libqt/qt.h"
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
     sprintf(space," ");
     fac = 0.0;
  }

  outfile->Printf("\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *                  %8s(T)                        *\n",name);
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf("\n");


  long int o = ndoccact;
  long int v = nvirt_no;

  long int oo   = o*o;
  long int vo   = v*o;
  long int vv   = v*v;
  long int ooo  = o*o*o;
  long int voo  = v*o*o;
  long int vvo  = v*v*o;
  long int vvv  = v*v*v;
  long int vooo = v*o*o*o;
  long int vvoo = v*v*o*o;

  double *F  = eps;
  double *E2ijak,**E2abci;
  E2ijak = (double*)malloc(vooo*sizeof(double));
  int nthreads = 1;
  #ifdef _OPENMP
      nthreads = Process::environment.get_n_threads();
  #endif

  long int memory = Process::environment.get_memory();
  if (options_["MEMORY"].has_changed()){
     memory  = options_.get_int("MEMORY");
     memory *= (long int)1024*1024;
  }
  // CDS // memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+3L*nthreads*v*v*v);
  long int memory_reqd = 8L*(2L*vvoo+vooo+vo+3L*nthreads*vvv);

  outfile->Printf("        num_threads:              %9i\n",nthreads);
  outfile->Printf("        available memory:      %9.2lf mb\n",(double)memory/1024./1024.);
  outfile->Printf("        memory requirements:   %9.2lf mb\n",
           (double)memory_reqd/1024./1024.);
  outfile->Printf("\n");


  long int nijk = 0;
  for (long int i=0; i<o; i++){
      for (long int j=0; j<=i; j++){
          for (long int k=0; k<=j; k++){
              nijk++;
          }
      }
  }
  long int**ijk = (long int**)malloc(nijk*sizeof(long int*));
  nijk = 0;
  for (long int i=0; i<o; i++){
      for (long int j=0; j<=i; j++){
          for (long int k=0; k<=j; k++){
              ijk[nijk] = (long int*)malloc(3*sizeof(long int));
              ijk[nijk][0] = i;
              ijk[nijk][1] = j;
              ijk[nijk][2] = k;
              nijk++;
          }
      }
  }
  outfile->Printf("        Number of ijk combinations: %ld\n",nijk);
  outfile->Printf("\n");


  E2abci = (double**)malloc(nthreads*sizeof(double*));
  // some v^3 intermediates
  double **Z  = (double**)malloc(nthreads*sizeof(double*));
  double **Z2 = (double**)malloc(nthreads*sizeof(double*));

  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(vvv*sizeof(double));
      Z[i]      = (double*)malloc(vvv*sizeof(double));
      Z2[i]     = (double*)malloc(vvv*sizeof(double));
  }

  std::shared_ptr<PSIO> psio(new PSIO());

  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&E2ijak[0],vooo*sizeof(double));
  psio->close(PSIF_DCC_IJAK,1);

  double *tempt = (double*)malloc(vvoo*sizeof(double));

  // first-order amplitudes for mp4
  if (ccmethod == 2) {
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],vvoo*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }

  for (long int a=0; a<vv; a++){
      C_DCOPY(oo,tb+a*oo,1,tempt+a,vv);
  }

  // might as well use t2's memory
  double*E2klcd = tb;
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb", (char*)&E2klcd[0],vvoo*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;
  outfile->Printf("        Computing (T) correction...\n");
  outfile->Printf("\n");
  outfile->Printf("        %% complete  total time\n");


  time_t stop,start = time(NULL);
  int pct10,pct20,pct30,pct40,pct50,pct60,pct70,pct80,pct90;
  pct10=pct20=pct30=pct40=pct50=pct60=pct70=pct80=pct90=0;

  /**
    *  if there is enough memory to explicitly thread, do so
    */
  #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
  for (long int ind=0; ind<nijk; ind++){
      long int i = ijk[ind][0];
      long int j = ijk[ind][1];
      long int k = ijk[ind][2];

      int thread = 0;
      #ifdef _OPENMP
          thread = omp_get_thread_num();
      #endif

      std::shared_ptr<PSIO> mypsio(new PSIO());
      mypsio->open(PSIF_DCC_ABCI,PSIO_OPEN_OLD);

      psio_address addr = psio_get_address(PSIO_ZERO,k*vvv*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],vvv*sizeof(double),addr,&addr);
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+j*vvo+i*vv,v,0.0,Z[thread],v*v);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+j*o*o*v+k*o*v,v,tempt+i*vvo,vv,1.0,Z[thread],v);

      //(ab)(ij)
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+i*vvo+j*vv,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+i*o*o*v+k*o*v,v,tempt+j*vvo,vv,1.0,Z2[thread],v);
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              C_DAXPY(v,1.0,Z2[thread]+b*vv+a*v,1,Z[thread]+a*vv+b*v,1);
          }
      }

      //(bc)(jk)
      addr = psio_get_address(PSIO_ZERO,(long int)j*vvv*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],vvv*sizeof(double),addr,&addr);
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+k*v*v*o+i*v*v,v,0.0,Z2[thread],v*v);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+k*voo+j*vo,v,tempt+i*vvo,vv,1.0,Z2[thread],v);
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              C_DAXPY(v,1.0,Z2[thread]+a*vv+b,v,Z[thread]+a*vv+b*v,1);
          }
      }

      //(ikj)(acb)
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+i*vvo+k*vv,v,0.0,Z2[thread],vv);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+i*voo+j*vo,v,tempt+k*vvo,vv,1.0,Z2[thread],v);
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              C_DAXPY(v,1.0,Z2[thread]+a*v+b,vv,Z[thread]+a*vv+b*v,1);
          }
      }

      //(ac)(ik)
      addr = psio_get_address(PSIO_ZERO,i*vvv*sizeof(double));
      mypsio->read(PSIF_DCC_ABCI,"E2abci",(char*)&E2abci[thread][0],vvv*sizeof(double),addr,&addr);
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+j*vvo+k*vv,v,0.0,Z2[thread],vv);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+j*voo+i*vo,v,tempt+k*vvo,vv,1.0,Z2[thread],v);
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              C_DAXPY(v,1.0,Z2[thread]+b*v+a,vv,Z[thread]+a*vv+b*v,1);
          }
      }

      //(ijk)(abc)
      F_DGEMM('t','t',vv,v,v,1.0,E2abci[thread],v,tempt+k*vvo+j*vv,v,0.0,Z2[thread],vv);
      F_DGEMM('n','t',v,vv,o,-1.0,E2ijak+k*voo+i*vo,v,tempt+j*vvo,vv,1.0,Z2[thread],v);
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              C_DAXPY(v,1.0,Z2[thread]+b*vv+a,v,Z[thread]+a*vv+b*v,1);
          }
      }

      C_DCOPY(vvv,Z[thread],1,Z2[thread],1);
      for (long int a=0; a<v; a++){
          double tai = t1[a*o+i];
          for (long int b=0; b<v; b++){
              long int ab = 1+(a==b);
              double tbj = t1[b*o+j];
              double E2iajb = E2klcd[i*vvo+a*vo+j*v+b];
              for (long int c=0; c<v; c++){
                  Z2[thread][a*vv+b*v+c] += fac*(tai * E2klcd[j*vvo+b*vo+k*v+c] +
                                                 tbj * E2klcd[i*vvo+a*vo+k*v+c] +
                                                 t1[c*o+k]*E2iajb);
                  Z2[thread][a*vv+b*v+c] /= (ab + (b==c) + (a==c));
              }
          }
      }

      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              for (long int c=0; c<v; c++){
                  long int abc = a*vv+b*v+c;
                  long int bac = b*vv+a*v+c;
                  long int acb = a*vv+c*v+b;
                  long int cba = c*vv+b*v+a;

                  E2abci[thread][abc] = Z2[thread][acb] + Z2[thread][bac] + Z2[thread][cba];
              }
          }
      }
      double dijk = F[i]+F[j]+F[k];
      long int ijkfac = ( 2-((i==j)+(j==k)+(i==k)) );
      // separate out these bits to save v^3 storage
      double tripval = 0.0;
      for (long int a=0; a<v; a++){
          double dijka = dijk-F[a+o];
          for (long int b=0; b<=a; b++){
              double dijkab = dijka-F[b+o];
              for (long int c=0; c<=b; c++){
                  long int abc = a*vv+b*v+c;
                  long int bca = b*vv+c*v+a;
                  long int cab = c*vv+a*v+b;
                  long int acb = a*vv+c*v+b;
                  long int bac = b*vv+a*v+c;
                  long int cba = c*vv+b*v+a;
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
      for (long int a=0; a<v; a++){
          for (long int b=0; b<v; b++){
              for (long int c=0; c<v; c++){
                  long int abc = a*vv+b*v+c;
                  long int bca = b*vv+c*v+a;
                  long int cab = c*vv+a*v+b;

                  E2abci[thread][abc]  = Z2[thread][abc] + Z2[thread][bca] + Z2[thread][cab];
              }
          }
      }
      tripval = 0.0;
      for (long int a=0; a<v; a++){
          double dijka = dijk-F[a+o];
          for (long int b=0; b<=a; b++){
              double dijkab = dijka-F[b+o];
              for (long int c=0; c<=b; c++){
                  long int abc = a*vv+b*v+c;
                  long int bca = b*vv+c*v+a;
                  long int cab = c*vv+a*v+b;
                  long int acb = a*vv+c*v+b;
                  long int bac = b*vv+a*v+c;
                  long int cba = c*vv+b*v+a;

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
            outfile->Printf("              %3.1lf  %8d s\n",100.0*ind/nijk,(int)stop-(int)start);

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
      outfile->Printf("\n");
      outfile->Printf("        (T) energy   %s                   %20.12lf\n",space,et);
      outfile->Printf("\n");
      outfile->Printf("        %s(T) correlation energy       %20.12lf\n",name,eccsd+et);
      outfile->Printf("      * %s(T) total energy             %20.12lf\n",name,eccsd+et+escf);
      outfile->Printf("\n");
  }else {
      emp4_t = myet;
      outfile->Printf("\n");
      outfile->Printf("        MP4(T) correlation energy:         %20.12lf\n",emp4_t);
      outfile->Printf("\n");
      outfile->Printf("        MP4(SDTQ) correlation energy:      %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+emp4_t);
      outfile->Printf("      * MP4(SDTQ) total energy:            %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+emp4_t+escf);
      outfile->Printf("\n");
  }


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
  delete[] name;
  delete[] space;

  return Success;
}

}} // end of namespaces

