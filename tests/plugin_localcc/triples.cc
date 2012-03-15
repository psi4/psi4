#include"ccsd.h"
#include"blas.h"
#ifdef _OPENMP
   #include<omp.h>
#endif

using namespace psi;

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

  double **tempRii = ccsd->wfn_->Rii->pointer();
  double *Rii = (double*)malloc(o*o*sizeof(double));
  for (int i=0; i<o; i++){
      for (int j=0; j<o; j++){
          Rii[i*o+j] = tempRii[i][j];
      }
  }
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
  if (options["CCMEMORY"].has_changed()){
     memory  = options.get_int("CCMEMORY");
     memory *= (long int)1024*1024;
  }
  memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+3L*nthreads*v*v*v);

  fprintf(outfile,"        num_threads =             %9i\n",nthreads);
  fprintf(outfile,"        available memory =     %9.2lf mb\n",memory/1024./1024.);
  fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
           8.*(2.*o*o*v*v+1.*o*o*o*v+(3.*nthreads)*v*v*v+1.*o*v)/1024./1024.);
  fprintf(outfile,"\n");
  fflush(outfile);

  bool threaded = true;
  if (memory<0){
     memory += (nthreads-1)*8L*3L*v*v*v;
     if (nthreads==1){
        fprintf(outfile,"        Error: not enough memory.\n");
        fprintf(outfile,"\n");
        fprintf(outfile,"        (T) requires at least %7.2lf mb\n",
             8.*(2.*o*o*v*v+1.*o*o*o*v+3.*v*v*v+1.*o*v)/1024./1024.);
        fprintf(outfile,"\n");
        fflush(outfile);
        return Failure;
     }
     threaded = false;
     nthreads = 1;
     fprintf(outfile,"        Not enough memory for explicit threading ... \n");
     fprintf(outfile,"\n");
     fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
              8.*(2.*o*o*v*v+1.*o*o*o*v+(3.)*v*v*v+1.*o*v)/1024./1024.);
     fprintf(outfile,"\n");
     fflush(outfile);
  }

  E2abci = (double**)malloc(nthreads*sizeof(double*));
  // some o^3 intermediates
  double **Z  = (double**)malloc(nthreads*sizeof(double*));
  double **Z2 = (double**)malloc(nthreads*sizeof(double*));
  double **Z3 = (double**)malloc(nthreads*sizeof(double*));
  double **Z4 = (double**)malloc(nthreads*sizeof(double*));

  boost::shared_ptr<PSIO> psio(new PSIO());
  double*tempE2=(double*)malloc(o*o*o*v*sizeof(double));
  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&tempE2[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);
  for (int i=0; i<o*o*o; i++){
      for (int a=0; a<v; a++){
          E2ijak[a*o*o*o+i] = tempE2[i*v+a];
      }
  }
  free(tempE2);

  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(v*v*v*sizeof(double));
      Z[i]      = (double*)malloc(o*o*o*sizeof(double));
      Z2[i]     = (double*)malloc(o*o*o*sizeof(double));
      Z3[i]     = (double*)malloc(o*o*o*sizeof(double));
      Z4[i]     = (double*)malloc(o*o*o*sizeof(double));
  }


  double *tempt = (double*)malloc(o*o*v*v*sizeof(double));

  if (ccsd->t2_on_disk){
     ccsd->tb = (double*)malloc(o*o*v*v*sizeof(double));
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&ccsd->tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }

  //for (int a=0; a<v*v; a++){
  //    F_DCOPY(o*o,ccsd->tb+a*o*o,1,tempt+a,v*v);
  //}
  F_DCOPY(o*o*v*v,ccsd->tb,1,tempt,1);

  // might as well use t2's memory
  double*E2klcd = ccsd->tb;
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd", (char*)&E2klcd[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;
  fprintf(outfile,"        Computing (T) correction...\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"        %% complete  total time\n");
  fflush(outfile);

  time_t stop,start = time(NULL);
  int pct10,pct20,pct30,pct40,pct50,pct60,pct70,pct80,pct90;
  pct10=pct20=pct30=pct40=pct50=pct60=pct70=pct80=pct90=0;

//heyheyhey
  int nabc = 0;
  for (int a=0; a<v; a++){
      for (int b=0; b<=a; b++){
          for (int c=0; c<=b; c++){
              nabc++;
          }
      }
  }
  int**abc = (int**)malloc(nabc*sizeof(int*));
  nabc = 0;
  for (int a=0; a<v; a++){
      for (int b=0; b<=a; b++){
          for (int c=0; c<=b; c++){
              abc[nabc] = (int*)malloc(3*sizeof(int));
              abc[nabc][0] = a;
              abc[nabc][1] = b;
              abc[nabc][2] = c;
              nabc++;
          }
      }
  }
  fprintf(outfile,"        Number of abc combinations: %i\n",nabc);
  fprintf(outfile,"\n");
  fflush(outfile);
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;

  /**
    *  if there is enough memory to explicitly thread, do so
    */
  if (threaded){
     #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
     for (int ind=0; ind<nabc; ind++){
         int a = abc[ind][0];
         int b = abc[ind][1];
         int c = abc[ind][2];

         int thread = 0;
         #ifdef _OPENMP
             thread = omp_get_thread_num();
         #endif

         boost::shared_ptr<PSIO> mypsio(new PSIO());
         mypsio->open(PSIF_ABCI4,PSIO_OPEN_OLD);
         psio_address addr = psio_get_address(PSIO_ZERO,(long int)(b*v*v*o+c*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
        
         // (1)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+a*o*o*v,o*o,0.0,Z[thread],o); 
         // (ikj)(acb)
         F_DGEMM('t','n',o,o*o,o,-1.0,tempt+c*o*o*v+a*o*o,o,E2ijak+b*o*o*o,o,1.0,Z[thread],o); 
         //for (int i=0; i<o; i++){
         //    for (int j=0; j<o; j++){
         //        for (int k=0; k<o; k++){
         //            double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][k*v+e] * tempt[a*o*o*v+e*o*o+i*o+j];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[b*o*o*o+i*o*o+j*o+m] * tempt[c*o*o*v+a*o*o+k*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += dum;
         //        }
         //    }
         //}

         addr = psio_get_address(PSIO_ZERO,(long int)(a*v*v*o+c*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ab)(ij)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+b*o*o*v,o*o,0.0,Z2[thread],o);
         //(ab)(ij)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+c*o*o*o,o,tempt+b*o*o*v+a*o*o,o,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 F_DAXPY(o,1.0,Z2[thread]+j*o*o+i*o,1,Z[thread]+i*o*o+j*o,1);
                 //for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][k*v+e] * tempt[b*o*o*v+e*o*o+j*o+i];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[c*o*o*o+i*o*o+k*o+m] * tempt[b*o*o*v+a*o*o+j*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+i*o+k];
                 //}
             }
         }

         addr = psio_get_address(PSIO_ZERO,(long int)(c*v*v*o+b*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(bc)(jk)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+a*o*o*v,o*o,0.0,Z2[thread],o);
         //(bc)(jk)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+b*o*o*o,o,tempt+a*o*o*v+c*o*o,o,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 F_DAXPY(o,1.0,Z2[thread]+i*o*o+j,o,Z[thread]+i*o*o+j*o,1);
                // for (int k=0; k<o; k++){
                //     double dum = 0.0;
                //     for (int e=0; e<v; e++){
                //         dum += E2abci[thread][j*v+e] * tempt[a*o*o*v+e*o*o+i*o+k];
                //     }
                //     for (int m=0; m<o; m++){
                //         dum -= E2ijak[b*o*o*o+k*o*o+j*o+m] * tempt[a*o*o*v+c*o*o+i*o+m];
                //     }
                //    Z[thread][i*o*o+j*o+k] += dum;
                // }
             }
         }
         addr = psio_get_address(PSIO_ZERO,(long int)(b*v*v*o+a*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ac)(ik)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+c*o*o*v,o*o,0.0,Z2[thread],o);
         //(ac)(ik)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+a*o*o*o,o,tempt+c*o*o*v+b*o*o,o,1.0,Z2[thread],o*o);
         //(1)
         F_DGEMM('t','t',o,o*o,o,-1.0,tempt+a*o*o*v+b*o*o,o,E2ijak+c*o*o*o,o*o,1.0,Z2[thread],o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][i*v+e] * tempt[c*o*o*v+e*o*o+k*o+j];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[a*o*o*o+j*o*o+i*o+m] * tempt[c*o*o*v+b*o*o+k*o+m];
                     //}
                     //  (1)
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[c*o*o*o+m*o*o+k*o+j] * tempt[a*o*o*v+b*o*o+i*o+m];
                     //}
                     Z[thread][i*o*o+j*o+k] += Z2[thread][k*o*o+j*o+i];
                 }
             }
         }
         addr = psio_get_address(PSIO_ZERO,(long int)(c*v*v*o+a*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ijk)(abc)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+b*o*o*v,o*o,0.0,Z2[thread],o);
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+a*o*o*o,o,tempt+b*o*o*v+c*o*o,o,1.0,Z2[thread],o*o);
         //(ijk)(abc)
         //for (int i=0; i<o; i++){
         //    for (int j=0; j<o; j++){
         //        for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][i*v+e] * tempt[b*o*o*v+e*o*o+j*o+k];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[a*o*o*o+k*o*o+i*o+m] * tempt[b*o*o*v+c*o*o+j*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+k*o+i];
         //        }
         //    }
         //}
         //(ikj)(acb)
         addr = psio_get_address(PSIO_ZERO,(long int)(a*v*v*o+b*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         F_DGEMM('n','n',o*o,o,v,1.0,tempt+c*o*o*v,o*o,E2abci[thread],v,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][j*v+e] * tempt[c*o*o*v+e*o*o+k*o+i];
                     //}
                     //Z[thread][i*o*o+j*o+k] += dum;
                     Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+k*o+i];
                 }
             }
         }

         
//heyheyhey
         // the singles part of (t) ...
         /*F_DCOPY(o*o*o,Z[thread],1,Z2[thread],1);
         for (int i=0; i<o; i++){
             double tai = t1[a*o+i];
             for (int j=0; j<o; j++){
                 double tbj = t1[b*o+j];
                 double E2iajb = E2klcd[i*v*v*o+a*v*o+j*v+b];
                 for (int k=0; k<o; k++){
                     Z2[thread][i*o*o+j*o+k] += (tai      *E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                                 tbj      *E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                                 t1[c*o+k]*E2iajb);
                 }
             }
         }
         double dabc = -F[a+o]-F[b+o]-F[c+o];
         for (int i=0; i<o; i++){
             double dabci = dabc + F[i];
             for (int j=0; j<o; j++){
                 double dabcij = dabci + F[j];
                 for (int k=0; k<o; k++){
                     double denom = dabcij + F[k];
                     Z[thread][i*o*o+j*o+k] /= denom;
                 }
             }
         }
         // transform Z and Z2 back to lmo basis
         F_DGEMM('n','t',o*o,o,o,1.0,Z[thread],o*o,Rii,o,0.0,E2abci[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     Z[thread][i*o*o+j*o+k] = 2.0*Z2[thread][i*o*o+j*o+k]
                                            - 3.0*Z2[thread][i*o*o+k*o+j]
                                            +     Z2[thread][j*o*o+k*o+i];
                     //E2abci[thread][i*o*o+j*o+k] = 2.0*Z2[thread][i*o*o+j*o+k]
                     //                       - 3.0*Z2[thread][i*o*o+k*o+j]
                     //                       +     Z2[thread][j*o*o+k*o+i];
                 }
             }
         }
         F_DGEMM('n','t',o*o,o,o,1.0,Z[thread],o*o,Rii,o,0.0,Z2[thread],o*o);

         // evaluate (t)
         for (int i=0; i<o; i++){
             double dum = 0.0;
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     dum +=  E2abci[thread][i*o*o+j*o+k] * Z2[thread][i*o*o+j*o+k];
                     //for (int ip=0; ip<o; ip++){
                     //for (int ipp=0; ipp<o; ipp++){
                     //    dum +=  E2abci[thread][i*o*o+j*o+ip] * Z[thread][i*o*o+j*o+ipp]*Rii[ip*o+k]*Rii[ipp*o+k]*ccsd->wfn_->centralfac[k];
                     //}
                     //}
                 }
             }
             etrip[thread] += dum * 2.0/3.0 * ccsd->wfn_->centralfac[i];
         }*/

//heyheyhey
         F_DCOPY(o*o*o,Z[thread],1,Z2[thread],1);
         double dabc = -F[a+o]-F[b+o]-F[c+o];
         for (int i=0; i<o; i++){
             double dabci = dabc+F[i];
             for (int j=0; j<o; j++){
                 double dabcij = dabci+F[j];
                 for (int k=0; k<o; k++){
                     double denom = dabcij+F[k];
                     Z[thread][i*o*o+j*o+k] /= denom;
                 }
             }
         }
         for (int i=0; i<o; i++){
             double tai = t1[a*o+i];
             for (int j=0; j<o; j++){
                 double tbj = t1[b*o+j];
                 double E2iajb = E2klcd[i*v*v*o+a*v*o+j*v+b];
                 for (int k=0; k<o; k++){
                     Z2[thread][i*o*o+j*o+k] += (tai      *E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                                 tbj      *E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                                 t1[c*o+k]*E2iajb);
                 }
             }
         }

         F_DCOPY(o*o*o,Z[thread],1,Z3[thread],1);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     Z3[thread][i*o*o+j*o+k] *= (1.0+0.5*(i==j)*(j==k));
                 }
             }
         }
         
         
         int abcfac = ( 2-((a==b)+(b==c)+(a==c)) );

         // transform one index of Z2 and Z3 to lmo basis:
         F_DGEMM('n','t',o*o,o,o,1.0,Z2[thread],o*o,Rii,o,0.0,Z4[thread],o*o);
         F_DCOPY(o*o*o,Z4[thread],1,E2abci[thread],1);
         F_DGEMM('n','t',o*o,o,o,1.0,Z3[thread],o*o,Rii,o,0.0,Z4[thread],o*o);

         // contribute to energy:
         double tripval = 0.0;
         for (int i=0; i<o; i++){
             double dum = 0.0;
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     dum         += Z4[thread][ijk] * E2abci[thread][ijk];
                     //dum         += Z2[thread][ijk] * Z3[thread][ijk];
                 }
             }
             tripval += dum * ccsd->wfn_->centralfac[i];
         }
         etrip[thread] += 3.0*tripval*abcfac;

         // Z3(ijk) = -2(Z(ijk) + jki + kij) + ikj + jik + kji
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int jki = j*o*o+k*o+i;
                     long int kij = k*o*o+i*o+j;
                     long int ikj = i*o*o+k*o+j;
                     long int jik = j*o*o+i*o+k;
                     long int kji = k*o*o+j*o+i;
                     Z3[thread][ijk] = -2.0*(Z[thread][ijk] + Z[thread][jki] + Z[thread][kij])
                                    +        Z[thread][ikj] + Z[thread][jik] + Z[thread][kji];
                 }
             }
         }

         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int ikj = i*o*o+k*o+j;
                     E2abci[thread][ijk] = Z2[thread][ikj]*0.5*(1.0+0.5*(i==j)*(j==k));
                 }
             }
         }

         // transform one index of of E2abci and Z3 to lmo basis:
         F_DGEMM('n','t',o*o,o,o,1.0,Z3[thread],o*o,Rii,o,0.0,Z4[thread],o*o);
         F_DCOPY(o*o*o,Z4[thread],1,Z3[thread],1);
         F_DGEMM('n','t',o*o,o,o,1.0,E2abci[thread],o*o,Rii,o,0.0,Z4[thread],o*o);

         // contribute to energy:
         tripval = 0.0;
         for (int i=0; i<o; i++){
             double dum = 0.0;
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     dum         += Z4[thread][ijk] * Z3[thread][ijk];
                     //dum         += E2abci[thread][ijk] * Z3[thread][ijk];
                 }
             }
             tripval += dum * ccsd->wfn_->centralfac[i];
         }
         etrip[thread] += tripval*abcfac;

         // the second bit
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     E2abci[thread][ijk] = Z2[thread][ijk]*0.5*(1.0+0.5*(i==j)*(j==k));
                 }
             }
         }

         // Z4 = Z(ijk)+jki+kij - 2( (ikj)+(jik)+(kji) )
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int jki = j*o*o+k*o+i;
                     long int kij = k*o*o+i*o+j;
                     long int ikj = i*o*o+k*o+j;
                     long int jik = j*o*o+i*o+k;
                     long int kji = k*o*o+j*o+i;
                     Z4[thread][ijk] =        Z[thread][ijk] + Z[thread][jki] + Z[thread][kij]
                                     - 2.0 * (Z[thread][ikj] + Z[thread][jik] + Z[thread][kji]);
                 }
             }
         }


         // transform one index of E2abci and Z4 to lmo basis:
         F_DGEMM('n','t',o*o,o,o,1.0,E2abci[thread],o*o,Rii,o,0.0,Z3[thread],o*o);
         F_DGEMM('n','t',o*o,o,o,1.0,Z4[thread],o*o,Rii,o,0.0,Z2[thread],o*o);

         // contribute to energy:
         tripval = 0.0;
         for (int i=0; i<o; i++){
             double dum = 0.0;
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     dum         += Z2[thread][ijk] * Z3[thread][ijk];
                     //dum         += E2abci[thread][ijk] * Z4[thread][ijk];
                 }
             }
             tripval += dum * ccsd->wfn_->centralfac[i];
         }
         etrip[thread] += tripval*abcfac;
//heyheyhey

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
               fprintf(outfile,"              %3.1lf  %8d s\n",100.0*ind/nabc,(int)stop-(int)start);
               fflush(outfile);
            }
         }
         mypsio->close(PSIF_ABCI4,1);
         mypsio.reset();
     }
  }
  else{
     fprintf(outfile,"on the to do pile!\n");
     return Failure;
  }


  double et = 0.0;
  for (int i=0; i<nthreads; i++) et += etrip[i];

  fprintf(outfile,"\n");
  if (ccsd->scale_t == 1.0){
     fprintf(outfile,"        (T) energy                   %20.12lf\n",et);
  }
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        (T) energy                   %20.12lf %20.12lf\n",et,et*ccsd->scale_t);
  }
  fprintf(outfile,"\n");
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf\n",ccsd->eccsd+et);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf %20.12lf\n",ccsd->eccsd+et,ccsd->eccsd+et*ccsd->scale_t);
  }
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"      * CCSD(T) total energy         %20.12lf\n",ccsd->eccsd+et+ccsd->escf);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"      * CCSD(T) total energy         %20.12lf %20.12lf\n",ccsd->eccsd+et+ccsd->escf,ccsd->eccsd+et*ccsd->scale_t+ccsd->escf);
  }
  fflush(outfile);
  ccsd->et = et;

  // free memory:
  free(Rii);
  if (ccsd->t2_on_disk){
     free(ccsd->tb);
  }
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
            
  return Success;
}


} // end of namespace



