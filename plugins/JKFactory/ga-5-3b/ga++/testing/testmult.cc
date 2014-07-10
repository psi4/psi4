#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ga++.h>

/* utilities for GA test programs */
#include "testutil.h"
#if HAVE_MATH_H
#   include <math.h>
#endif

#define N 256          /* first dimension  */
#define BASE 0
#define PERMUTE_ 

#define GA_DATA_TYPE MT_C_FLOAT
#define GA_ABS(a) (((a) >= 0) ? (a) : (-(a)))
#define TOLERANCE 0.000001

#ifdef MPIPP
#include <mpi.h>
#define CLOCK_ MPI_Wtime
#else
#include "tcgmsg.h"
#define CLOCK_ tcg_time
#endif

DoublePrecision gTime=0.0, gStart;

void
test(int data_type, int ndim) {
  int me=GA_Nodeid();
  GA::GlobalArray *g_a, *g_b, *g_c, *g_A, *g_B, *g_C;
  int dims[GA_MAX_DIM]={N,N,2,2,2,1,1};
  int lo[GA_MAX_DIM]={1,1,1,1,1,0,0};
  int hi[GA_MAX_DIM]={N-2,N-2,1,1,1,0,0};
  int clo[2], chi[2], m, n, k;
  double value1_dbl = 2.0, value2_dbl = 2.0;
  double alpha_dbl = 1.0, beta_dbl = 0.0;
  float value1_flt = 2.0, value2_flt = 2.0;
  float alpha_flt = 1.0, beta_flt = 0.0;
  DoubleComplex value1_dcpl = {2.0, 2.0}, value2_dcpl = {2.0, 2.0};
  DoubleComplex alpha_dcpl = {1.0, 0.0} , beta_dcpl = {0.0, 0.0}; 
  void *value1, *value2, *alpha, *beta;
  char name_a[] = "array A";
  char name_b[] = "array B";
  char name_c[] = "array C";
  char name_a_[] = "array A_";
  char name_b_[] = "array B_";
  char name_c_[] = "array C_";
  
  switch (data_type) {
  case C_FLOAT:
    alpha  = (void *)&alpha_flt;
    beta   = (void *)&beta_flt;
    value1 = (void *)&value1_flt;
    value2 = (void *)&value2_flt;
    if(me==0) printf("Single Precision: Testing GA_Sgemm,NGA_Matmul_patch for %d-Dimension", ndim);
    break;      
  case C_DBL:
    alpha  = (void *)&alpha_dbl;
    beta   = (void *)&beta_dbl;
    value1 = (void *)&value1_dbl;
    value2 = (void *)&value2_dbl;
    if(me==0) printf("Double Precision: Testing GA_Dgemm,NGA_Matmul_patch for %d-Dimension", ndim); 
    break;    
  case C_DCPL:
    alpha  = (void *)&alpha_dcpl;
    beta   = (void *)&beta_dcpl;
    value1 = (void *)&value1_dcpl;
    value2 = (void *)&value2_dcpl;
    if(me==0) printf("Double Complex:   Testing GA_Zgemm,NGA_Matmul_patch for %d-Dimension", ndim);
    break;
  default:
    GA::SERVICES.error("wrong data type", data_type);
  }

  g_a = GA::SERVICES.createGA(data_type, ndim, dims, name_a, NULL);
  g_b = GA::SERVICES.createGA(g_a, name_b);  
  g_c = GA::SERVICES.createGA(g_a, name_c);

  g_a->fill(value1);
  g_b->fill( value2);
  g_c->zero();
  
  /** g_c = g_a * g_b */
  g_c->matmulPatch('N', 'N', alpha, beta,
		   g_a, lo, hi,
		   g_b, lo, hi,
		   lo, hi);  
  g_a->destroy();
  g_b->destroy();
  
  /** 
   * Verifying g_c:
   * 1. Create g_A(=g_a) and g_B(=g_b)
   * 2. g_C = g_A*g_B; (Using Gemm routines)
   * 3. g_A = g_c; (copy the 2-d patch og g_c into g_A)
   * 4. g_C = g_A - g_C; (Using add() routine by making beta=-1.0)
   * 5. If all the elements in g_C is zero, implies SUCCESS.
   */
  dims[0] = dims[1] = m = n = k = N-2; 
  g_A = GA::SERVICES.createGA(data_type, 2, dims, name_a_, NULL);
  g_B = GA::SERVICES.createGA(g_A, name_b_);
  g_C = GA::SERVICES.createGA(g_A, name_c_);

  g_A->fill(value1);
  g_B->fill(value2);
  g_C->zero();

  gStart = CLOCK_();
  switch (data_type) {
  case C_FLOAT:
    g_C->sgemm('N', 'N', m, n, k, alpha_flt,  g_A, g_B, beta_flt);
    beta_flt = -1.0;
    break;      
  case C_DBL:
    g_C->dgemm('N', 'N', m, n, k, alpha_dbl,  g_A, g_B, beta_dbl);
    beta_dbl = -1.0;
    break;    
  case C_DCPL:
    g_C->zgemm('N', 'N', m, n, k, alpha_dcpl, g_A, g_B, beta_dcpl);
    beta_dcpl.real = -1.0; 
    break;
  default:
    GA::SERVICES.error("wrong data type", data_type);
  }
  gTime += CLOCK_()-gStart;

  g_B->destroy();
  
  clo[0] = clo[1] = 0;
  chi[0] = chi[1] = N-3;

  g_A->copyPatch('N', g_c, lo, hi, clo, chi) ;  
  g_C->add(alpha, g_A, beta, g_C);

  switch (data_type) {
  case C_FLOAT:
    value1_flt = g_C->fdot(g_C);
    if(value1_flt != 0.0)
      GA::SERVICES.error("GA_Sgemm, NGA_Matmul_patch Failed", 0);
    break;
  case C_DBL:
    value1_dbl = g_C->ddot(g_C);
    if(value1_dbl != 0.0)
      GA::SERVICES.error("GA_Dgemm, NGA_Matmul_patch Failed", 0);
    break;
  case C_DCPL:
    value1_dcpl = g_C->zdot(g_C);
    if(value1_dcpl.real != 0.0 || value1_dcpl.imag != 0.0)
      GA::SERVICES.error("GA_Zgemm, NGA_Matmul_patch Failed", 0);
    break;
  default:
    GA::SERVICES.error("wrong data type", data_type);
  }  
  
  if(me==0) printf("....OK\n");

  g_A->destroy();
  g_c->destroy();
  g_C->destroy();
}

void
do_work() {
  int i;
  int me = GA_Nodeid();
 
  for(i=2; i<=GA_MAX_DIM; i++) {
     test(C_FLOAT, i);
     test(C_DBL,   i);
     test(C_DCPL,  i);
     if(me == 0) printf("\n\n");
    GA::SERVICES.sync();
  }
}
     

int 
main(int argc, char **argv) {

Integer heap=9000000, stack=9000000;
int me;
DoublePrecision time;

 GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
 me=GA_Nodeid();

 time = CLOCK_(); 
 do_work();

#ifdef TIME
 printf("%d Total Time = %lf\n", me, CLOCK_()-time);
 printf("%d GEMM Total Time = %lf\n", me, gTime);
#endif
 
 if(me==0)printf("\nSuccessfull\n\n");
 GA::Terminate();
 return 0;
}

