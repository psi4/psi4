/**
 * @file testmatmultc.c
 * @author Jeff Daily, PNNL jeff.daily@pnl.gov
 *
 * This is the C version of testmatmul.F.  The port was as direct as possible,
 * hopefully.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <math.h>
#include <stdlib.h>

#include "mp3.h"
#include "ga.h"
#include "macdecls.h"
#include "xgemm.h"


void load_ga(int handle, double *f, int dim1, int dim2);
void verify_ga_dgemm(char xt1, char xt2, int num_m, int num_n, int num_k,
        double alpha, int g_a, int g_b, double beta, int g_c,
        double *tmpa, double *tmpb, double *tmpc);


#define dgemm_verify 1
#define nummax 1024
#define howmany 2
#define ntrans 4


/*
 * test ga_dgemm
 * Note: - change nummax for large arrays
 *       - turn off "dgemm_verify" for large arrays due to memory 
 *         limitations, as dgemm_verify=1 for large arrays produces 
 *         segfault, dumps core,or any crap.
 */
int main(int argc, char **argv)
{
    int num_m;
    int num_n;
    int num_k;
    int i;
    int ii;
    double *h0;
    int g_c;
    int g_b;
    int g_a;
    double a;
    double t1;
    double mf;
    double avg_t[ntrans];
    double avg_mf[ntrans];
    int itime;
    int ntimes;
    int me;
    int nums_m[/*howmany*/] = {512,1024};
    int nums_n[/*howmany*/] = {512,1024};
    int nums_k[/*howmany*/] = {512,1024};
    char transa[/*ntrans*/] = "ntnt";
    char transb[/*ntrans*/] = "nntt";
    char ta;
    char tb;
    double *tmpa;
    double *tmpb;
    double *tmpc;
    int ndim;
    int dims[2];
#ifdef BLOCK_CYCLIC
    int block_size[2];
#endif

    MP_INIT(argc,argv);
    if (!MA_init(MT_DBL,1,20000000)) {
        GA_Error("failed: ma_init(MT_DBL,1,20000000)",10);
    }
    GA_INIT(argc,argv);
    me = GA_Nodeid();

    h0 = (double*)malloc(sizeof(double) * nummax*nummax);
    tmpa = (double*)malloc(sizeof(double) * nummax*nummax);
    tmpb = (double*)malloc(sizeof(double) * nummax*nummax);
    tmpc = (double*)malloc(sizeof(double) * nummax*nummax);

    ii = 0;
    for (i=0; i<nummax*nummax; i++) {
        ii = ii + 1;
        if (ii > nummax) {
            ii = 0;
        }
        h0[i] = ii;
    }

    /* Compute times assuming 500 mflops and 5 second target time */
    /* ntimes = max(3.0d0,5.0d0/(4.0d-9*num**3)); */
    ntimes = 5;

    for (ii=0; ii<howmany; ii++) {
        num_m = nums_m[ii];
        num_n = nums_n[ii];
        num_k = nums_k[ii];
        a = 0.5/(num_m*num_n);
        if (num_m > nummax || num_n > nummax || num_k > nummax) {
            GA_Error("Insufficient memory: check nummax", 1);
        }

#ifndef BLOCK_CYCLIC
        ndim = 2;

        dims[0] = num_m;
        dims[1] = num_n;
        if (!((g_c = NGA_Create(MT_DBL,ndim,dims,"g_c",NULL)))) {
            GA_Error("failed: create g_c",20);
        }

        dims[0] = num_k;
        dims[1] = num_n;
        if (!((g_b = NGA_Create(MT_DBL,ndim,dims,"g_b",NULL)))) {
            GA_Error("failed: create g_b",30);
        }

        dims[0] = num_m;
        dims[1] = num_k;
        if (!((g_a = NGA_Create(MT_DBL,ndim,dims,"g_a",NULL)))) {
            GA_Error("failed: create g_a",40);
        }
#else
        ndim = 2;
        block_size[0] = 128;
        block_size[1] = 128;

        dims[0] = num_m;
        dims[1] = num_n;
        g_c = GA_Create_handle();
        GA_Set_data(g_c,ndim,dims,MT_DBL);
        GA_Set_array_name(g_c,"g_c");
        GA_Set_block_cyclic(g_c,block_size);
        if (!GA_Allocate(g_c)) {
            GA_Error("failed: create g_c",40);
        }

        dims[0] = num_k;
        dims[1] = num_n;
        g_b = GA_Create_handle();
        GA_Set_data(g_b,ndim,dims,MT_DBL);
        GA_Set_array_name(g_b,"g_b");
        GA_Set_block_cyclic(g_b,block_size);
        if (!ga_allocate(g_b)) {
            GA_Error("failed: create g_b",40);
        }

        dims[0] = num_m;
        dims[1] = num_k;
        g_a = GA_Create_handle();
        GA_Set_data(g_a,ndim,dims,MT_DBL);
        GA_Set_array_name(g_a,"g_a");
        GA_Set_block_cyclic(g_a,block_size);
        if (!ga_allocate(g_a)) {
            GA_Error('failed: create g_a',40);
        }
#endif         

        /* Initialize matrices A and B */
        if (me == 0) { 
            load_ga(g_a, h0, num_m, num_k);
            load_ga(g_b, h0, num_k, num_n);
        }
        GA_Zero(g_c);
        GA_Sync();

        if (GA_Nodeid() == 0) {
            printf("\nMatrix Multiplication on C = A[%ld,%ld]xB[%ld,%ld]\n",
                    (long)num_m, (long)num_k, (long)num_k, (long)num_n);
            fflush(stdout);
        }

        for (i=0; i<ntrans; i++) {
            avg_t[i]  = 0.0;
            avg_mf[i] = 0.0;
        }

        for (itime=0; itime<ntimes; itime++) {
            for (i=0; i<ntrans; i++) {
                GA_Sync();
                ta = transa[i];
                tb = transb[i];
                t1 = MP_TIMER();
                GA_Dgemm(ta,tb,num_m,num_n,num_k,1.0, g_a, g_b, 0.0, g_c);
                t1 = MP_TIMER() - t1;
                if (me == 0) {
                    mf = 2e0*num_m*num_n*num_k/t1*1e-6/GA_Nnodes();
                    avg_t[i]  = avg_t[i]+t1;
                    avg_mf[i] = avg_mf[i] + mf;
                    printf("%15s%2d: %12.4f seconds %12.1f mflops/proc  %c %c\n",
                            "Run#", itime, t1, mf, ta, tb);
                    fflush(stdout);
                    if (dgemm_verify && itime == 0) {
                        /* recall the C API swaps the matrix order */
                        /* we swap it here for the Fortran-based verify */
                        verify_ga_dgemm(tb, ta, num_n, num_m, num_k, 1.0,
                                g_b, g_a, 0.0, g_c, tmpb, tmpa, tmpc);
                    }
                }
            }
        }

        if (GA_Nodeid() == 0) {
            printf("\n");
            for (i=0; i<ntrans; i++) {
                printf("%17s: %12.4f seconds %12.1f mflops/proc  %c %c\n",
                        "Average", avg_t[i]/ntimes, avg_mf[i]/ntimes,
                        transa[i], transb[i]);
            }
            if(dgemm_verify) {
                printf("All GA_Dgemms are verified...O.K.\n");
            }
            fflush(stdout);
        }

        /*
           GA_Print(g_a);
           GA_Print(g_b);
           GA_Print(g_c);
           */

        GA_Destroy(g_c);
        GA_Destroy(g_b);
        GA_Destroy(g_a);
    }

    /* ???
       format(a15, i2, ': ', e12.4, ' seconds ',f12.1, 
       .     ' mflops/proc ', 3a2)
       */
    if (GA_Nodeid() == 0) {
        printf("All tests successful\n");
    }

    free(h0);
    free(tmpa);
    free(tmpb);
    free(tmpc);

    GA_Terminate();
    MP_FINALIZE();
    return 0;
}


/*
 * Verify for correctness. Process 0 computes BLAS dgemm 
 * locally. For larger arrays, disbale this test as memory
 * might not be sufficient
 */
void verify_ga_dgemm(char xt1, char xt2, int num_m, int num_n, int num_k,
        double alpha, int g_a, int g_b, double beta, int g_c,
        double *tmpa, double *tmpb, double *tmpc)
{
    int i,j,type,ndim,dims[2],lo[2],hi[2];
    double abs_value;

    for (i=0; i<num_n; i++) {
        for (j=0; j<num_m; j++) {
            tmpc[j+i*num_m] = -1.0;
            tmpa[j+i*num_m] = -2.0;
        }
    }

    NGA_Inquire(g_a, &type, &ndim, dims);
    lo[0] = 0;
    lo[1] = 0;
    hi[0] = dims[0]-1;
    hi[1] = dims[1]-1;
    NGA_Get(g_a, lo, hi, tmpa, &dims[1]);

    NGA_Inquire(g_a, &type, &ndim, dims);
    lo[0] = 0;
    lo[1] = 0;
    hi[0] = dims[0]-1;
    hi[1] = dims[1]-1;
    NGA_Get(g_b, lo, hi, tmpb, &dims[1]);

    /* compute dgemm sequentially */
    xb_dgemm(&xt1, &xt2, &num_m, &num_n, &num_k,
            &alpha, tmpa, &num_m,
            tmpb, &num_k, &beta,
            tmpc, &num_m);

    /* after computing c locally, verify it with the values in g_c */
    NGA_Inquire(g_a, &type, &ndim, dims);
    lo[0] = 0;
    lo[1] = 0;
    hi[0] = dims[0]-1;
    hi[1] = dims[1]-1;
    NGA_Get(g_c, lo, hi, tmpa, &dims[1]);

    for (i=0; i<num_n; i++) {
        for (j=0; j<num_m; j++) {
            abs_value = fabs(tmpc[j+i*num_m]-tmpa[j+i*num_m]);
            if(abs_value > 1.0 || abs_value < -1.0) {
                printf("Values are = %f %f\n",
                        tmpc[j+i*num_m], tmpa[j+i*num_m]);
                printf("Values are = %f %f\n", 
                        fabs(tmpc[j+i*num_m]-tmpa[j*i*num_m]), abs_value);
                fflush(stdout);
                GA_Error("verify ga_dgemm failed", 1);
            }
        }
    }
}

/**
 * called by process '0' (or your master process )
 */
void load_ga(int handle, double *f, int dim1, int dim2)
{
      int lo[2], hi[2];
      
      if (dim1 < 0 || dim2 < 0) {
          return;
      }

      lo[0] = 0;
      lo[1] = 0;
      hi[0] = dim1-1;
      hi[1] = dim2-1;
      NGA_Put(handle, lo, hi, f, &dim1);
}

/*
c
c-----------------------------------------------------------------------
c     must be called by all processors, if you need to fillup the 
c     entire array
c
      subroutine load_ga_from_square(handle,num,f,ndim)
      implicit none
      integer handle, memhandle
      integer num,ndim
      real*8 f(ndim,ndim)
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ga_nodeid, i1, i2, i, j, ix, jx
#include "mafdecls.fh"

      call ga_distribution(handle, ga_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

c     nx = ihi - ilo + 1
c     ny = jhi - jlo + 1

      do i = ilo,ihi,ndim
         do j = jlo,jhi,ndim
            call ga_put(handle,i,min(ihi,i+ndim),j,min(jhi,j+ndim),
     &                    f,ndim)
         enddo
      enddo

      return
      end
*/

/*
c
c-----------------------------------------------------------------------
c     must be called by all processors, if you need to fillup the 
c     entire array
c
      subroutine load_ga_from_triangle(handle,f,ndim)
      implicit none
      integer handle, memhandle
      real*8 f(*)
      integer ndim
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ga_nodeid, i1, i2, i, j, ix, jx
#include "mafdecls.fh"

      call ga_distribution(handle, ga_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

      nx = ihi - ilo + 1
      ny = jhi - jlo + 1

      if (.not.ma_alloc_get(MT_DBL,nx*ny,'flap',memhandle,ibuff)) then
         call ga_error('failed: allocate triangle',100)
      endif

      do i = 1,nx
         do j = 1,ny
            ix = i + ilo - 1
            jx = j + jlo - 1
            i1 = min(ix,jx)
            i2 = max(ix,jx)
            dbl_mb(ibuff + nx*(j-1) + (i-1) ) = f(i2*(i2-1)/2 + i1)
         enddo
      enddo

      call ga_put(handle,ilo,ihi,jlo,jhi,
     &              dbl_mb(ibuff),nx)

      if (.not.ma_free_heap(memhandle)) then
         call ga_error('failed: free triangle',100)
      endif

      return
      end
*/
