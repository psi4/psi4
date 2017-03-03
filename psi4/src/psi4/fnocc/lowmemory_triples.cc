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

PsiReturnType CoupledCluster::lowmemory_triples() {

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


  outfile->Printf("\n");
  outfile->Printf( "        Warning: due to limited available memory,\n");
  outfile->Printf( "        using less efficient, low-memory algorithm\n");
  outfile->Printf("\n");

  long int o = ndoccact;
  long int v = nvirt_no;

  long int oo   = o*o;
  long int vo   = v*o;
  long int ooo  = o*o*o;
  long int voo  = v*o*o;
  long int vvo  = v*v*o;
  long int vooo = v*o*o*o;
  long int vvoo = v*v*o*o;

  double *F  = eps;
  double *E2ijak,**E2abci;
  // CDS // E2ijak = (double*)malloc(o*o*o*v*sizeof(double));
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
  // CDS // memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+5L*nthreads*o*o*o);
  long int memory_reqd = 8L*(2L*vvoo+vooo+vo+5L*nthreads*ooo);

  outfile->Printf("        num_threads:              %9i\n",nthreads);
  outfile->Printf("        available memory:      %9.2lf mb\n",(double)memory/1024./1024.);
  outfile->Printf("        memory requirements:   %9.2lf mb\n",
           (double)memory_reqd/1024./1024.);
  outfile->Printf("\n");


  bool threaded = true;

  /*
  if (memory_reqd > memory){
     // CDS // memory += (nthreads-1)*8L*5L*ooo;
     if (nthreads==1){
        outfile->Printf("        Error: not enough memory.\n");
        outfile->Printf("\n");
        outfile->Printf("        (T) requires at least %7.2lf mb\n",
             (double)(2.*o*o*v*v+1.*o*o*o*v+5.*o*o*o+1.*o*v)/1024./1024.);
        outfile->Printf("\n");

        return Failure;
     }
     threaded = false;
     nthreads = 1;
     outfile->Printf("        Not enough memory for explicit threading ... \n");
     outfile->Printf("\n");
     outfile->Printf("        memory requirements =  %9.2lf mb\n",
              8.*(2.*o*o*v*v+1.*o*o*o*v+(5.)*o*o*o+1.*o*v)/1024./1024.);
     outfile->Printf("\n");

  }
  */

  // CDS updated
  if (memory_reqd > memory) {
     outfile->Printf("        Not enough memory for requested threading ...\n");
     outfile->Printf("\n");

     long int min_memory_reqd = 8L*(2L*vvoo+vooo+vo+5L*ooo);
     if (min_memory_reqd > memory) {
         outfile->Printf("        Sorry, not even enough memory for 1 thread.\n");
	 delete [] name;
	 delete [] space;
	 free(E2ijak);
         return Failure;
     }

     long int mem_leftover = memory - min_memory_reqd;
     int extra_threads = (int) (mem_leftover / 5L*ooo);
     nthreads = 1 + extra_threads;
     outfile->Printf("        Attempting to proceed with %d threads\n",
       nthreads);
  }

  E2abci = (double**)malloc(nthreads*sizeof(double*));
  // some o^3 intermediates
  double **Z  = (double**)malloc(nthreads*sizeof(double*));
  double **Z2 = (double**)malloc(nthreads*sizeof(double*));
  double **Z3 = (double**)malloc(nthreads*sizeof(double*));
  double **Z4 = (double**)malloc(nthreads*sizeof(double*));

  std::shared_ptr<PSIO> psio(new PSIO());
  double*tempE2=(double*)malloc(vooo*sizeof(double));
  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&tempE2[0],vooo*sizeof(double));
  psio->close(PSIF_DCC_IJAK,1);
  for (long int i=0; i<ooo; i++){
      for (long int a=0; a<v; a++){
          E2ijak[a*ooo+i] = tempE2[i*v+a];
      }
  }
  free(tempE2);

  long int dim = ooo > vo ? ooo : vo;
  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(dim*sizeof(double));
      Z[i]      = (double*)malloc(ooo*sizeof(double));
      Z2[i]     = (double*)malloc(ooo*sizeof(double));
      Z3[i]     = (double*)malloc(ooo*sizeof(double));
      Z4[i]     = (double*)malloc(ooo*sizeof(double));
  }

  double *tempt = (double*)malloc(vvoo*sizeof(double));

  if (t2_on_disk){
     // CDS // tb = (double*)malloc((long int)o*o*v*v*sizeof(double));
     tb = (double*)malloc(vvoo*sizeof(double));
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tb[0],vvoo*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }
  if (ccmethod == 2) {
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],vvoo*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }

  C_DCOPY(vvoo,tb,1,tempt,1);

  // might as well use t2's memory
  double*E2klcd = tb;
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb", (char*)&E2klcd[0],
      vvoo*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;

  time_t stop,start = time(NULL);
  int pct10,pct20,pct30,pct40,pct50,pct60,pct70,pct80,pct90;
  pct10=pct20=pct30=pct40=pct50=pct60=pct70=pct80=pct90=0;

  long int nabc = 0;
  for (long int a=0; a<v; a++){
      for (long int b=0; b<=a; b++){
          for (long int c=0; c<=b; c++){
              nabc++;
          }
      }
  }
  long int**abc = (long int**)malloc(nabc*sizeof(long int*));
  nabc = 0;
  for (long int a=0; a<v; a++){
      for (long int b=0; b<=a; b++){
          for (long int c=0; c<=b; c++){
              abc[nabc] = (long int*)malloc(3*sizeof(long int));
              abc[nabc][0] = a;
              abc[nabc][1] = b;
              abc[nabc][2] = c;
              nabc++;
          }
      }
  }
  outfile->Printf("        Number of abc combinations: %i\n",nabc);
  outfile->Printf("\n");

  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;

  outfile->Printf("        Computing (T) correction...\n");
  outfile->Printf("\n");
  outfile->Printf("        %% complete  total time\n");

  /**
    *  if there is enough memory to explicitly thread, do so
    */

  std::vector< std::shared_ptr<PSIO> > mypsio;
  for (int i = 0; i < nthreads; i++) {
      mypsio.push_back( (std::shared_ptr<PSIO>)(new PSIO()) );
      mypsio[i]->open(PSIF_DCC_ABCI4,PSIO_OPEN_OLD);
  }

  if (threaded){
     #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
     for (long int ind=0; ind<nabc; ind++){
         long int a = abc[ind][0];
         long int b = abc[ind][1];
         long int c = abc[ind][2];

         int thread = 0;
         #ifdef _OPENMP
             thread = omp_get_thread_num();
         #endif

         //std::shared_ptr<PSIO> mypsio(new PSIO());
         //mypsio->open(PSIF_DCC_ABCI4,PSIO_OPEN_OLD);
         psio_address addr = psio_get_address(PSIO_ZERO,(b*vvo+c*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],vo*sizeof(double),addr,&addr);

         // (1)
         F_DGEMM('t','t',o,oo,v,1.0,E2abci[thread],v,tempt+a*voo,oo,0.0,Z[thread],o);
         // (ikj)(acb)
         F_DGEMM('t','n',o,oo,o,-1.0,tempt+c*voo+a*oo,o,E2ijak+b*ooo,o,1.0,Z[thread],o);

         addr = psio_get_address(PSIO_ZERO,(a*vvo+c*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],vo*sizeof(double),addr,&addr);
         //(ab)(ij)
         F_DGEMM('t','t',o,oo,v,1.0,E2abci[thread],v,tempt+b*voo,oo,0.0,Z2[thread],o);
         //(ab)(ij)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+c*ooo,o,tempt+b*voo+a*oo,o,1.0,Z2[thread],oo);
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 C_DAXPY(o,1.0,Z2[thread]+j*oo+i*o,1,Z[thread]+i*oo+j*o,1);
             }
         }

         addr = psio_get_address(PSIO_ZERO,(c*vvo+b*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],vo*sizeof(double),addr,&addr);
         //(bc)(jk)
         F_DGEMM('t','t',o,oo,v,1.0,E2abci[thread],v,tempt+a*voo,oo,0.0,Z2[thread],o);
         //(bc)(jk)
         F_DGEMM('t','n',oo,o,o,-1.0,E2ijak+b*ooo,o,tempt+a*voo+c*oo,o,1.0,Z2[thread],oo);
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 C_DAXPY(o,1.0,Z2[thread]+i*oo+j,o,Z[thread]+i*oo+j*o,1);
             }
         }
         addr = psio_get_address(PSIO_ZERO,(b*vvo+a*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],vo*sizeof(double),addr,&addr);
         //(ac)(ik)
         F_DGEMM('t','t',o,oo,v,1.0,E2abci[thread],v,tempt+c*voo,oo,0.0,Z2[thread],o);
         //(ac)(ik)
         F_DGEMM('t','n',oo,o,o,-1.0,E2ijak+a*ooo,o,tempt+c*voo+b*oo,o,1.0,Z2[thread],oo);
         //(1)
         F_DGEMM('t','t',o,oo,o,-1.0,tempt+a*voo+b*oo,o,E2ijak+c*ooo,oo,1.0,Z2[thread],o);
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     Z[thread][i*oo+j*o+k] += Z2[thread][k*oo+j*o+i];
                 }
             }
         }
         addr = psio_get_address(PSIO_ZERO,(c*vvo+a*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],vo*sizeof(double),addr,&addr);
         //(ijk)(abc)
         F_DGEMM('t','t',o,oo,v,1.0,E2abci[thread],v,tempt+b*voo,oo,0.0,Z2[thread],o);
         F_DGEMM('t','n',oo,o,o,-1.0,E2ijak+a*ooo,o,tempt+b*voo+c*oo,o,1.0,Z2[thread],oo);
         //(ijk)(abc)
         //(ikj)(acb)
         addr = psio_get_address(PSIO_ZERO,(a*vvo+b*vo)*sizeof(double));
         mypsio[thread]->read(PSIF_DCC_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         F_DGEMM('n','n',oo,o,v,1.0,tempt+c*voo,oo,E2abci[thread],v,1.0,Z2[thread],oo);
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     Z[thread][i*oo+j*o+k] += Z2[thread][j*oo+k*o+i];
                 }
             }
         }

         C_DCOPY(ooo,Z[thread],1,Z2[thread],1);
         double dabc = -F[a+o]-F[b+o]-F[c+o];
         for (long int i=0; i<o; i++){
             double dabci = dabc+F[i];
             for (long int j=0; j<o; j++){
                 double dabcij = dabci+F[j];
                 for (long int k=0; k<o; k++){
                     double denom = dabcij+F[k];
                     Z[thread][i*oo+j*o+k] /= denom;
                 }
             }
         }
         for (long int i=0; i<o; i++){
             double tai = t1[a*o+i];
             for (long int j=0; j<o; j++){
                 double tbj = t1[b*o+j];
                 double E2iajb = E2klcd[i*vvo+a*vo+j*v+b];
                 for (long int k=0; k<o; k++){
                     Z2[thread][i*oo+j*o+k] += fac *
                         (tai * E2klcd[j*vvo+b*vo+k*v+c] +
                          tbj * E2klcd[i*vvo+a*vo+k*v+c] +
                          t1[c*o+k]*E2iajb);
                 }
             }
         }

         C_DCOPY(ooo,Z[thread],1,Z3[thread],1);
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     Z3[thread][i*oo+j*o+k] *= (1.0+0.5*(i==j)*(j==k));
                 }
             }
         }


         long int abcfac = ( 2-((a==b)+(b==c)+(a==c)) );

         // contribute to energy:
         double tripval = 0.0;
         for (long int i=0; i<o; i++){
             double dum = 0.0;
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     dum         += Z3[thread][ijk] * Z2[thread][ijk];
                 }
             }
             tripval += dum;
         }
         etrip[thread] += 3.0*tripval*abcfac;

         // Z3(ijk) = -2(Z(ijk) + jki + kij) + ikj + jik + kji
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     long int jki = j*oo+k*o+i;
                     long int kij = k*oo+i*o+j;
                     long int ikj = i*oo+k*o+j;
                     long int jik = j*oo+i*o+k;
                     long int kji = k*oo+j*o+i;
                     Z3[thread][ijk] = -2.0*(Z[thread][ijk] + Z[thread][jki] + Z[thread][kij])
                                    +        Z[thread][ikj] + Z[thread][jik] + Z[thread][kji];
                 }
             }
         }

         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     long int ikj = i*oo+k*o+j;
                     E2abci[thread][ijk] = Z2[thread][ikj]*0.5*(1.0+0.5*(i==j)*(j==k));
                 }
             }
         }

         // contribute to energy:
         tripval = 0.0;
         for (long int i=0; i<o; i++){
             double dum = 0.0;
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     dum         += E2abci[thread][ijk] * Z3[thread][ijk];
                 }
             }
             tripval += dum;
         }
         etrip[thread] += tripval*abcfac;

         // the second bit
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     E2abci[thread][ijk] = Z2[thread][ijk]*0.5*(1.0+0.5*(i==j)*(j==k));
                 }
             }
         }

         // Z4 = Z(ijk)+jki+kij - 2( (ikj)+(jik)+(kji) )
         for (long int i=0; i<o; i++){
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     long int jki = j*oo+k*o+i;
                     long int kij = k*oo+i*o+j;
                     long int ikj = i*oo+k*o+j;
                     long int jik = j*oo+i*o+k;
                     long int kji = k*oo+j*o+i;
                     Z4[thread][ijk] =        Z[thread][ijk] + Z[thread][jki] + Z[thread][kij]
                                     - 2.0 * (Z[thread][ikj] + Z[thread][jik] + Z[thread][kji]);
                 }
             }
         }

         // contribute to energy:
         tripval = 0.0;
         for (long int i=0; i<o; i++){
             double dum = 0.0;
             for (long int j=0; j<o; j++){
                 for (long int k=0; k<o; k++){
                     long int ijk = i*oo+j*o+k;
                     dum         += Z4[thread][ijk] * E2abci[thread][ijk];
                 }
             }
             tripval += dum;
         }
         etrip[thread] += tripval*abcfac;

         // print out update
         if (thread==0){
            int print = 0;
            stop = time(NULL);
            if ((double)ind/nabc >= 0.1 && !pct10){      pct10 = 1; print=1;}
            else if ((double)ind/nabc >= 0.2 && !pct20){ pct20 = 1; print=1;}
            else if ((double)ind/nabc >= 0.3 && !pct30){ pct30 = 1; print=1;}
            else if ((double)ind/nabc >= 0.4 && !pct40){ pct40 = 1; print=1;}
            else if ((double)ind/nabc >= 0.5 && !pct50){ pct50 = 1; print=1;}
            else if ((double)ind/nabc >= 0.6 && !pct60){ pct60 = 1; print=1;}
            else if ((double)ind/nabc >= 0.7 && !pct70){ pct70 = 1; print=1;}
            else if ((double)ind/nabc >= 0.8 && !pct80){ pct80 = 1; print=1;}
            else if ((double)ind/nabc >= 0.9 && !pct90){ pct90 = 1; print=1;}
            if (print){
               outfile->Printf("              %3.1lf  %8d s\n",100.0*ind/nabc,(int)stop-(int)start);

            }
         }
         //mypsio->close(PSIF_DCC_ABCI4,1);
         //mypsio.reset();
     }
  }
  else{
     outfile->Printf("on the to do pile!\n");
     delete [] name;
     delete [] space;
     free(E2ijak);
     for (int i=0; i<nthreads; i++) free(Z4[i]);
     free(Z4);
     free(etrip);
     nabc = 0;
     for (long int a=0; a<v; a++)
       for (long int b=0; b<=a; b++)
	 for (long int c=0; c<=b; c++)
	   free(abc[nabc++]);
     free(abc);

     return Failure;
  }
  for (int i = 0; i < nthreads; i++) {
      mypsio[i]->close(PSIF_DCC_ABCI4,1);
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


  delete[] name;
  delete[] space;

  // free memory:
  free(E2ijak);
  for (int i=0; i<nthreads; i++){
      free(E2abci[i]);
      free(Z[i]);
      free(Z2[i]);
      free(Z3[i]);
      free(Z4[i]);
  }
  free(Z);
  free(Z2);
  free(Z3);
  free(Z4);
  free(E2abci);
  free(etrip);

  nabc = 0;
  for (long int a=0; a<v; a++)
    for (long int b=0; b<=a; b++)
      for (long int c=0; c<=b; c++)
	free(abc[nabc++]);
  free(abc);


  return Success;
}


}} // end of namespaces

