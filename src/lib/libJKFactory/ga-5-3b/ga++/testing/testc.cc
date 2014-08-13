#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;
#include "ga++.h"

#define N 5
#define GA_DATA_TYPE MT_F_REAL

/*
using std::cout;
using std::printf;
using std::sin;
using std::endl;
*/

int
main(int argc, char *argv[]) {
  
  int ONE=1;   /* useful constants */
  int n=N, type=MT_F_DBL;
  int me, nproc;
  int i, row;
  int dims[2]={N,N};
  int lo[2], hi[2];
  
  /* Note: on all current platforms DoublePrecision == double */
  double buf[N], err, alpha, beta;

  int heap  = 200000;
  int stack = 200000;
  
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  me = GA_Nodeid();
  nproc = GA_Nnodes();
  
  if(me==0)printf("Size: %d: Creating matrix A\n", nproc);
  GA::GlobalArray *g_a = GA::SERVICES.createGA(type, 2, dims, (char *)"A", NULL);
 
  
  if(me==0)printf("OK\n");
  
  if(me==0)printf("Creating matrix B\n");
  /* create matrix B  so that it has dims and distribution of A*/
  GA::GlobalArray *g_b = GA::SERVICES.createGA(g_a, (char *)"B");
  if(me==0)printf("OK\n");
  
  g_a->zero();   /* zero the matrix */
  
  if(me==0)printf("Initializing matrix A\n");
  /* fill in matrix A with random values in range 0.. 1 */ 
  lo[1]=0; hi[1]=n-1;
  for(row=me; row<n; row+= nproc){
    /* each process works on a different row in MIMD style */
    lo[0]=hi[0]=row;   
    for(i=0; i<n; i++) buf[i]=sin((double)i + 0.1*(row+1));
    g_a->put(lo, hi, buf, &n);
  }
//  g_a->print();
  
  
  if(me==0)printf("Symmetrizing matrix A\n");
  g_a->symmetrize();   /* symmetrize the matrix A = 0.5*(A+A') */
  
  /* check if A is symmetric */ 
  if(me==0)printf("Checking if matrix A is symmetric\n");
  g_a->transpose(g_b); /* B=A' */
  alpha=1.; beta=-1.;
  g_b->add(&alpha, g_a, &beta, g_b);  /* B= A - B */
  err= g_b->ddot(g_b);
  
  if(me==0)printf("Error=%f\n",(double)err);
  
  if(me==0)printf("\nChecking atomic accumulate \n");
  	
  g_a->zero();   /* zero the matrix */
  for(i=0; i<n; i++) buf[i]=(double)i;
  
  /* everybody accumulates to the same location/row */
  alpha = 1.0;
  row = n/2;
  lo[0]=hi[0]=row;
  lo[1]=0; hi[1]=n-1;
  g_a->acc(lo, hi, buf, &ONE, &alpha );
  GA::SERVICES.sync();
  
  if(me==0){ /* node 0 is checking the result */
    
    g_a->get(lo, hi, buf,&ONE);
    for(i=0; i<n; i++) if(buf[i] != (double)nproc*i)
      GA::SERVICES.error((char *)"failed: column=",i);
    printf("OK\n\n");
    
  }

  g_a->destroy();
  g_b->destroy();

  if(me==0) cout << "Terminating...\n";
  GA::Terminate();
  
}
