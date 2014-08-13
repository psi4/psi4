#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "macdecls.h"
#include "ga.h"
#include "mp3.h"

/* utilities for GA test programs */
#include "testutil.h"

#define N 400          /* first dimension  */
#define BASE 0
#define PERMUTE_ 

#define GA_DATA_TYPE MT_C_FLOAT
#define GA_ABS(a) (((a) >= 0) ? (a) : (-(a)))
#define TOLERANCE 0.000001


DoublePrecision gTime=0.0, gStart;

void
test(int data_type, int ndim) {
  int me=GA_Nodeid();
  int g_a, g_b, g_c, g_A, g_B, g_C;
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
  SingleComplex value1_scpl = {2.0, 2.0}, value2_scpl = {2.0, 2.0};
  SingleComplex alpha_scpl = {1.0, 0.0} , beta_scpl = {0.0, 0.0}; 
  void *value1=NULL, *value2=NULL, *alpha=NULL, *beta=NULL;

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
  case C_SCPL:
    alpha  = (void *)&alpha_scpl;
    beta   = (void *)&beta_scpl;
    value1 = (void *)&value1_scpl;
    value2 = (void *)&value2_scpl;
    if(me==0) printf("Single Complex:   Testing GA_Cgemm,NGA_Matmul_patch for %d-Dimension", ndim);
    break;
  default:
    GA_Error("wrong data type", data_type);
  }

  g_a = NGA_Create(data_type, ndim, dims, "array A", NULL);
  g_b = GA_Duplicate(g_a, "array B");  
  g_c = GA_Duplicate(g_a, "array C");
  if(!g_a || !g_b || !g_c) GA_Error("Create failed: a, b or c",1);

  GA_Fill(g_a, value1);
  GA_Fill(g_b, value2);
  GA_Zero(g_c);

  NGA_Matmul_patch('N', 'N', alpha, beta,
           g_a, lo, hi,
           g_b, lo, hi,
           g_c, lo, hi);  
  GA_Destroy(g_a);
  GA_Destroy(g_b);

  /** 
   * Verifying g_c:
   * 1. Create g_A(=g_a) and g_B(=g_b)
   * 2. g_C = g_A*g_B; (Using Gemm routines)
   * 3. g_A = g_c; (copy the 2-d patch og g_c into g_A)
   * 4. g_C = g_A - g_C; (Using add() routine by making beta=-1.0)
   * 5. If all the elements in g_C is zero, implies SUCCESS.
   */
  dims[0] = dims[1] = m = n = k = N-2; 
  g_A = NGA_Create(data_type, 2, dims, "array A_", NULL);
  g_B = GA_Duplicate(g_A, "array B_");
  g_C = GA_Duplicate(g_A, "array C_");
  if(!g_A || !g_B || !g_C) GA_Error("Create failed: A, B or C",n);
  GA_Fill(g_A, value1);
  GA_Fill(g_B, value2);
  GA_Zero(g_C);

  gStart = MP_TIMER();
  switch (data_type) {
  case C_FLOAT:
    GA_Sgemm('N', 'N', m, n, k, alpha_flt,  g_A, g_B, beta_flt, g_C);
    beta_flt = -1.0;
    break;      
  case C_DBL:
    GA_Dgemm('N', 'N', m, n, k, alpha_dbl,  g_A, g_B, beta_dbl, g_C);
    beta_dbl = -1.0;
    break;    
  case C_DCPL:
    GA_Zgemm('N', 'N', m, n, k, alpha_dcpl, g_A, g_B, beta_dcpl, g_C);
    beta_dcpl.real = -1.0; 
    break;
  case C_SCPL:
    GA_Cgemm('N', 'N', m, n, k, alpha_scpl, g_A, g_B, beta_scpl, g_C);
    beta_scpl.real = -1.0; 
    break;
  default:
    GA_Error("wrong data type", data_type);
  }
  gTime += MP_TIMER()-gStart;

  GA_Destroy(g_B);
  
  clo[0] = clo[1] = 0;
  chi[0] = chi[1] = N-3;

  NGA_Copy_patch('N', g_c, lo, hi, g_A, clo, chi) ;
  
  GA_Add(alpha, g_A, beta, g_C, g_C);
  /*  NGA_Add_patch (alpha, g_c, lo, hi, beta, g_C, clo, chi, g_C, clo, chi);*/

  switch (data_type) {
  case C_FLOAT:
    value1_flt = GA_Fdot(g_C, g_C);
    if(fabsf(value1_flt) > TOLERANCE) {
      printf("\nabs(result) = %f > %f\n", fabsf(value1_flt), TOLERANCE);
      GA_Error("GA_Sgemm, NGA_Matmul_patch Failed", 1);
    }
    break;
  case C_DBL:
    value1_dbl = GA_Ddot(g_C, g_C);
    if(fabs(value1_dbl) > TOLERANCE) {
      printf("\nabs(result) = %f > %f\n", fabs(value1_dbl), TOLERANCE);
      GA_Error("GA_Dgemm, NGA_Matmul_patch Failed", 1);
    }
    break;
  case C_DCPL:
    value1_dcpl = GA_Zdot(g_C, g_C);
    if(fabs(value1_dcpl.real) > TOLERANCE
            || fabs(value1_dcpl.imag) > TOLERANCE) {
      printf("\nabs(result) = %f+%fi > %f\n",
              fabs(value1_dcpl.real), fabs(value1_dcpl.imag), TOLERANCE);
      GA_Error("GA_Zgemm, NGA_Matmul_patch Failed", 1);
    }
    break;
  case C_SCPL:
    value1_scpl = GA_Cdot(g_C, g_C);
    if(fabsf(value1_scpl.real) > TOLERANCE
            || fabsf(value1_scpl.imag) > TOLERANCE) {
      printf("\nabs(result) = %f+%fi > %f\n",
              fabsf(value1_scpl.real), fabsf(value1_scpl.imag), TOLERANCE);
      GA_Error("GA_Sgemm, NGA_Matmul_patch Failed", 1);
    }
    break;
  default:
    GA_Error("wrong data type", data_type);
  }  
  
  if(me==0) printf("....OK\n");

  GA_Destroy(g_A);
  GA_Destroy(g_c);
  GA_Destroy(g_C);
}

void
do_work() {
  int i;
  int me = GA_Nodeid();
 
  for(i=2; i<=GA_MAX_DIM; i++) {
     test(C_FLOAT, i);
     test(C_DBL,   i);
     test(C_DCPL,  i);
     test(C_SCPL,  i);
     if(me == 0) printf("\n\n");
    GA_Sync();
  }
}
     

int 
main(int argc, char **argv) {

Integer heap=9000000, stack=9000000;
int me, nproc;
DoublePrecision time;

    MP_INIT(argc,argv);

    GA_INIT(argc,argv);                           /* initialize GA */

    nproc = GA_Nnodes();
    me = GA_Nodeid();

    if(me==0) printf("Using %d processes\n\n",nproc);

    if(!MA_init((Integer)MT_F_DBL, stack/nproc, heap/nproc))
       GA_Error("MA_init failed bytes= %d",stack+heap);   

#ifdef PERMUTE
      {
        int i, *list = (int*)malloc(nproc*sizeof(int));
        if(!list)GA_Error("malloc failed",nproc);

        for(i=0; i<nproc;i++)list[i]=nproc-1-i;

        GA_Register_proclist(list, nproc);
        free(list);
      }
#endif

    if(GA_Uses_fapi())GA_Error("Program runs with C API only",1);

    time = MP_TIMER();
    do_work();
    /*    printf("%d: Total Time = %lf\n", me, MP_TIMER()-time);
      printf("%d: GEMM Total Time = %lf\n", me, gTime);
    */

    if(me==0)printf("\nSuccess\n\n");
    GA_Terminate();

    MP_FINALIZE();

    return 0;
}

