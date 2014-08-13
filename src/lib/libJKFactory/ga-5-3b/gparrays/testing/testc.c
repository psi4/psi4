/* Test algorithm parameters */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <math.h>
#include "ga.h"
#include "gp.h"
#include "macdecls.h"
#include "mp3.h"

#include <stdlib.h>

/*
#define N_I  4
#define N_J  4
*/
#define N_I  32
#define N_J  32
#define N_K  32
/*
#define Q_I 2
#define Q_J 2
*/
#define Q_I 8
#define Q_J 8

/* get random patch for GP */
void get_range( int ndim, int dims[], int lo[], int hi[], int g_p)
{
  int dim, nproc, i, itmp;
  nproc = NGA_Nodeid();
  /* Mix up values on different processors */
  for (i=0; i<nproc; i++) itmp = rand();

  for(dim=0; dim<ndim; dim++){
    int toss1, toss2;
    toss1 = rand()%dims[dim];
    toss2 = rand()%dims[dim];
    if (toss1<toss2) {
      lo[dim]=toss1;
      hi[dim]=toss2;
    } else {
      hi[dim]=toss1;
      lo[dim]=toss2;
    }
  }
  /*
  nproc = (nproc+8)%NGA_Nnodes();
  GP_Distribution(g_p, nproc, lo, hi);
  */
}

void factor(int p,int idim, int jdim, int kdim, int *pdi, int *pdj, int *pdk)
{
  int MAX_FACTOR = 10000;
  int ii, i, j, k, ip, ifac, pmax, ichk;
  int ti, tj, tk;
  int *prime, *fac;

  prime = malloc(MAX_FACTOR*sizeof(int));
  fac = malloc(MAX_FACTOR*sizeof(int));
  /**
   * factor p completely.
   * First find all prime numbers, besides 1, less than
   * or equal to sqrt of p
   */
  ip = (int)sqrt((double)p) + 1;
  pmax = 0;
  for (i=2; i<=ip; i++) {
    ichk = 1;
    for (j=0; j<pmax; j++) {
      if (i%prime[j] == 0) {
        ichk = 0;
        break;
      }
    }
    if (ichk == 1) {
      pmax++;
      if (pmax > MAX_FACTOR) printf("Overflow in grid factor\n");
      prime[pmax-1] = i;
    }
  }
  
  /**
   *  find all prime factors of p
   */
  ip = p;
  ifac = 0;
  for (i=0; i<pmax; i++) {
    while (ip%prime[i] == 0) {
      fac[ifac] = prime[i];
      ifac++;
      ip = ip/prime[i];
    }
  }

  ti = 1;
  tj = 1;
  tk = 1;
  for (ii=ifac-1; ii>=0; ii--) {
    i = idim/ti;
    j = jdim/tj;
    k = kdim/tk;
    if (i >= j && i >= k && i > 1) {
      ti = fac[ii]*ti;
    } else if (j >= i && j >= k && j > 1) {
      tj = fac[ii]*tj;
    } else if (k >= i && k >= j && k > 1) {
      tk = fac[ii]*tk;
    } else {
      printf("Too many processors in factoring routine\n");
    }
  }

  free(prime);
  free(fac);

  *pdi = ti;
  *pdj = tj;
  *pdk = tk;
}

void do_work()
{
  int g_p, me, i, ii, j, jj, l, k, nv;
  int nproc, next;
  int m_k_ij, m_l_ij, idx;
  int dims[2], lo[2], hi[2], ndim;
  int lo_t[2], hi_t[2];
  int dims3[3], lo3[3], hi3[3], blocks[3], chunk[3];
  int nelems, nsize;
  int idim, jdim, kdim, subscript[2], size;
  int pdi, pdj, pdk;
  int ld[2], ld_sz[2];
  int checksize;
  int *ptr, *mapc;
  void **buf_ptr;
  void *buf;
  int *buf_size;
  void *elem_buf;
  int *subscripts;

  /* Create Global Pointer array */
  dims[0] = N_I;
  dims[1] = N_J;
  ndim = 2;
  me = NGA_Nodeid();
  nproc = NGA_Nnodes();

  g_p = GP_Create_handle();
  GP_Set_dimensions(g_p, ndim, dims);
  GP_Allocate(g_p);

  /* Find locally owned elements in Global Pointer array.
     Only these elements can be assigned to data using the
     GP_Assign_local_element routine. */
  GP_Distribution(g_p, me, lo, hi);
  idim = hi[0] - lo[0] + 1;
  jdim = hi[1] - lo[1] + 1;
  for (ii=0; ii<idim; ii++) {
    i = ii + lo[0];
    for (jj=0; jj<jdim; jj++) {
      j = jj + lo[1];
      idx = j*N_I + i;
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      /* Allocate local memory for object and assign it values */
      size = sizeof(int)*(m_k_ij*m_l_ij+2);
      ptr = (int*)GP_Malloc(size);
      ptr[0] = m_k_ij;
      ptr[1] = m_l_ij;
      for (k=0; k<m_k_ij; k++) {
        for (l=0; l<m_l_ij; l++) {
          ptr[l*m_k_ij+k+2] = l*m_k_ij+k+idx;
        }
      }
      subscript[0] = i;
      subscript[1] = j;
      
      if (subscript[0]<lo[0] || subscript[0]>hi[0] || subscript[1]<lo[1] ||
          subscript[1]>hi[1]) {
        printf("p[%d] assign i: %d j: %d lo[0]: %d hi[0]: %d lo[1]: %d hi[1]: %d\n",
            me,subscript[0],subscript[1],lo[0],hi[0],lo[1],hi[1]);
      }
      GP_Assign_local_element(g_p, subscript, (void*)ptr, size);
    }
  }
  
  /* Guarantee data consistency */
  GP_Sync();
  /*GP_Debug(g_p);*/

  /* Generate bounding coordinates to an arbitrary patch in GP array */
  get_range(ndim, dims, lo, hi, g_p);

  /* Find the total amount of data contained in the patch */
  nsize = (hi[0]-lo[0]+1)*(hi[1]-lo[1]+1);
  GP_Get_size(g_p, lo, hi, &size);

  /* Allocate local buffers and retrieve data */
  buf = (void*)malloc(size);
  buf_ptr = (void**)malloc(nsize*sizeof(void*));
  buf_size = (int*) malloc(nsize*sizeof(int));
  ld[1] = hi[0]-lo[0]+1;
  ld[0] = hi[1]-lo[1]+1;
  ld_sz[1] = hi[0]-lo[0]+1;
  ld_sz[0] = hi[1]-lo[1]+1;
  GA_Set_debug(1);
  GP_Get(g_p, lo, hi, buf, buf_ptr, ld, buf_size, ld_sz, &size, 0);
  if (me==0) printf("\nCompleted GP_Get\n");
  GA_Set_debug(0);
  
  /* Check contents of buffers to see if data is as expected */
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      idx = j*N_I + i;
      ptr = (int*)buf_ptr[ii*ld[0]+jj];
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      if (buf_size[ii*ld_sz[0]+jj] != 4*(ptr[0]*ptr[1]+2)) {
        printf("p[%d] size expected: %d actual: %d\n",me,buf_size[ii*ld_sz[0]+jj],
            4*(ptr[0]*ptr[1]+2));
      }
      if (ptr[0] != m_k_ij) {
        printf("p[%d] [%d,%d] Dimension(1) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
      }
      if (ptr[1] != m_l_ij) {
        printf("p[%d] [%d,%d] Dimension(1) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
      }
      for (k=0; k<ptr[0]; k++) {
        for (l=0; l<ptr[1]; l++) {
          if (ptr[l*ptr[0]+k+2] != l*m_k_ij+k+idx) {
            printf("p[%d] Element(1) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
                me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
          }
        }
      }
    }
  }
  if (me==0) printf("\nCompleted check of GP_Get\n");

  /* Clear local buffers */
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      ptr = (int*)buf_ptr[ii*ld[0]+jj];
      size = ptr[0]*ptr[1]+2;
      for (k=0; k<size; k++) {
        ptr[k] = 0;
      }
    }
  }

  /* Get data using information on buffers */
  GP_Get(g_p, lo, hi, buf, buf_ptr, ld, buf_size, ld_sz, &size, 1);

  /* Recheck contents of buffers to see if data is as expected */
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      idx = j*N_I + i;
      ptr = (int*)buf_ptr[ii*ld[0]+jj];
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      if (buf_size[ii*ld_sz[0]+jj] != 4*(ptr[0]*ptr[1]+2)) {
        printf("p[%d] size expected: %d actual: %d\n",me,buf_size[ii*ld_sz[0]+jj],
            4*(ptr[0]*ptr[1]+2));
      }
      if (ptr[0] != m_k_ij) {
        printf("p[%d] [%d,%d] Dimension(1) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
      }
      if (ptr[1] != m_l_ij) {
        printf("p[%d] [%d,%d] Dimension(1) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
      }
      for (k=0; k<ptr[0]; k++) {
        for (l=0; l<ptr[1]; l++) {
          if (ptr[l*ptr[0]+k+2] != l*m_k_ij+k+idx) {
            printf("p[%d] Element(1) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
                me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
          }
        }
      }
    }
  }
  if (me==0) printf("\nCompleted check of GP_Get using known buffer sizes\n");
  GP_Sync();

  /* Clear all bits in GP_Array */
  GP_Memzero(g_p);
  /* Test to see if all bits actually are zero */
  GP_Distribution(g_p, me, lo_t, hi_t);
  for (i=lo_t[0]; i<=hi_t[0]; i++) {
    ii = i - lo_t[0];
    subscript[0] = i;
    for (j=lo_t[1]; j<=hi_t[1]; j++) {
      jj = j - lo_t[1];
      subscript[1] = j;
      GP_Access_element(g_p, subscript, &elem_buf, &size);
      ptr = (int*)elem_buf;
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      if (size/4 != m_k_ij*m_l_ij+2) {
        printf("p[%d] Mismatched sizes in memzero test\n",me);
      }
      for (k=0; k<m_k_ij*m_l_ij+2; k++) {
        if (ptr[k] != 0) {
          printf("p[%d] Nonzero element %d in memzero test ptr[%d]: %c\n",
              me,idx,k,ptr[k]);
        }
      }
    }
  }
  GP_Sync();
  if (me==0) printf("\nZeroed all bits in GP array\n");

  /* dellocate buffers and resize them to test put capability */
  free(buf);
  free(buf_ptr);
  free(buf_size);

  next = (me+1)%nproc;
  GP_Distribution(g_p, me, lo, hi);
  GP_Get_size(g_p, lo, hi, &size);

  nelems = (hi[0]-lo[0]+1)*(hi[1]-lo[1]+1);
  ld[0] = hi[1]-lo[1]+1;
  ld_sz[0] = hi[1]-lo[1]+1;
  buf = (void*)malloc(size);
  buf_ptr = (void**)malloc(nelems*sizeof(void*));
  buf_size = (int*) malloc(nelems*sizeof(int));

  /* Fill buffers with contents of GP array */
  nsize = 0;
  next = 0;
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      idx = j*N_I + i;
      ptr = (int*)(((char*)buf)+next);
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      nsize = sizeof(int)*(m_k_ij*m_l_ij+2);
      next += nsize;
      buf_size[ii*ld_sz[0] +jj] = nsize;
      buf_ptr[ii*ld[0]+jj] = ptr;
      ptr[0] = m_k_ij;
      ptr[1] = m_l_ij;
      for (k=0; k<m_k_ij; k++) {
        for (l=0; l<m_l_ij; l++) {
          ptr[l*m_k_ij + k + 2] = l*m_k_ij + k + idx;
        }
      }
    }
  }
  checksize = 1;
  GP_Put(g_p, lo, hi, buf_ptr, ld, buf_size, ld_sz, &size, checksize);
  GP_Sync();

  /* Check contents of GP array */
  GP_Distribution(g_p, me, lo, hi);
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    subscript[0] = i;
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      subscript[1] = j;
      idx = j*N_I + i;
      GP_Access_element(g_p, subscript, &elem_buf, &size);
      ptr = (int*)elem_buf;
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      if (size != sizeof(int)*(m_k_ij*m_l_ij+2)) {
        printf("p[%d] [%d,%d] Size actual: %d expected: %ld\n",
               me,i,j,size,(long)(sizeof(int)*(m_k_ij*m_l_ij+2)));
      }
      if (ptr[0] != m_k_ij) {
        printf("p[%d] [%d,%d] Dimension(2) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
      }
      if (ptr[1] != m_l_ij) {
        printf("p[%d] [%d,%d] Dimension(2) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
      }
      for (k=0; k<ptr[0]; k++) {
        for (l=0; l<ptr[1]; l++) {
          if (ptr[l*m_k_ij+k+2] != l*m_k_ij+k+idx) {
            printf("p[%d] Element(2) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
                me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
          }
        }
      }
    }
  }
  if (me==0) printf("\nCompleted check of GP_Put\n");

  /* Test gather. Deallocate buffers first */
  free(buf);
  free(buf_ptr);
  free(buf_size);

  nv = 5;
  subscripts = (int*)malloc(2*nv*sizeof(int));
  for (ii=0; ii<nv; ii++) {
    jj = me + ii*(N_I-1);
    i = jj%N_I;
    j = (jj-i)/N_I;
    subscripts[ii*2] = i;
    subscripts[ii*2+1] = j;
  }

  /* Get size of data to be gathered */
  GP_Gather_size(g_p, nv, subscripts, &size);

  /* Check size for correct value */
  idx = 0;
  for (ii=0; ii<nv; ii++) {
    i = subscripts[ii*2];
    j = subscripts[ii*2+1];
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    idx += sizeof(int)*(m_k_ij*m_l_ij+2);
  }
  if (idx != size) {
    printf("p[%d] Size actual: %d expected: %d\n",me,size,idx);
  }
  if (me == 0) printf("\nCompleted check of GP_Gather_size\n");

  buf = (void*)malloc(size);
  buf_ptr = (void**)malloc(nv*sizeof(void*));
  buf_size = (int*) malloc(nv*sizeof(int));

  /* Gather data elements */
  GP_Gather(g_p, nv, subscripts, buf, buf_ptr, buf_size, &size, 0);

  /* Check data in buffers to see if it is correct */
  for (ii=0; ii<nv; ii++) {
    i = subscripts[ii*2];
    j = subscripts[ii*2+1];
    idx = j*N_I + i;
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    ptr = (int*)buf_ptr[ii];
    if ((int)buf_size[ii] != sizeof(int)*(m_k_ij*m_l_ij+2)) {
      printf("p[%d] [%d,%d] Size(3) i actual: %d expected: %ld\n",
          me,i,j,buf_size[ii],(long)(sizeof(int)*(m_k_ij*m_l_ij+2)));
    }
    if (ptr[0] != m_k_ij) {
      printf("p[%d] [%d,%d] Dimension(3) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
    }
    if (ptr[1] != m_l_ij) {
      printf("p[%d] [%d,%d] Dimension(3) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
    }
    for (k=0; k<ptr[0]; k++) {
      for (l=0; l<ptr[1]; l++) {
        if (ptr[l*m_k_ij+k+2] != l*m_k_ij+k+idx) {
          printf("p[%d] Element(3) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
              me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
        }
      }
    }
  }
  if (me==0) printf("\nCompleted check of GP_Gather\n");

  /* clean contents of buffers */
  for (ii=0; ii<nv; ii++) {
    i = subscripts[ii*2];
    j = subscripts[ii*2+1];
    idx = j*N_I + i;
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    ptr = (int*)buf_ptr[ii];
    size = ptr[0]*ptr[1]+2;
    for (k=0; k<size; k++){
      ptr[k] = 0;
    }
  }

  GP_Gather(g_p, nv, subscripts, buf, buf_ptr, buf_size, &size, 1);

  /* Check data in buffers to see if it is correct */
  for (ii=0; ii<nv; ii++) {
    i = subscripts[ii*2];
    j = subscripts[ii*2+1];
    idx = j*N_I + i;
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    ptr = (int*)buf_ptr[ii];
    if ((int)buf_size[ii] != sizeof(int)*(m_k_ij*m_l_ij+2)) {
      printf("p[%d] [%d,%d] Size(3) i actual: %d expected: %ld\n",
          me,i,j,buf_size[ii],(long)(sizeof(int)*(m_k_ij*m_l_ij+2)));
    }
    if (ptr[0] != m_k_ij) {
      printf("p[%d] [%d,%d] Dimension(3) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
    }
    if (ptr[1] != m_l_ij) {
      printf("p[%d] [%d,%d] Dimension(3) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
    }
    for (k=0; k<ptr[0]; k++) {
      for (l=0; l<ptr[1]; l++) {
        if (ptr[l*m_k_ij+k+2] != l*m_k_ij+k+idx) {
          printf("p[%d] Element(3) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
              me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
        }
      }
    }
  }
  free(subscripts);
  if (me==0) printf("\nCompleted check of GP_Gather using known buffer sizes\n");
  

  /* Clean up buffers and clear all bits in GP_Array */
  GP_Sync();
  free(buf);
  free(buf_ptr);
  free(buf_size);
  GP_Memzero(g_p);
  GP_Sync();

  /* Test scatter capability. Start by figuring out how many elements processor
   * will scatter and how large the total data size is
   */
  size = 0;
  nsize = dims[0]*dims[1];
  nelems = 0;
  for (ii = me; ii < nsize; ii += nproc) {
    nelems++;
    i = ii%dims[0];
    j = (ii-i)/dims[0];
    idx = j*N_I + i;
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    size += sizeof(int)*(m_k_ij*m_l_ij+2);
  }
  buf = (void*)malloc(size);
  buf_ptr = (void**)malloc(nelems*sizeof(void*));
  buf_size = (int*)malloc(nelems*sizeof(int));
  subscripts = (int*)malloc(2*nelems*sizeof(int));

  /* Set up buf_ptr and buf_size arrays */
  elem_buf = buf;
  nv = 0;
  for (ii = me; ii < nsize; ii += nproc) {
    i = ii%dims[0];
    j = (ii-i)/dims[0];
    idx = j*N_I + i;
    m_k_ij = i%Q_I + 1;
    m_l_ij = j%Q_J + 1;
    size = sizeof(int)*(m_k_ij*m_l_ij+2);
    subscripts[nv*2] = i;
    subscripts[nv*2+1] = j;
    buf_ptr[nv] = elem_buf;
    buf_size[nv] = size;
    ptr = (int*)buf_ptr[nv];
    ptr[0] = m_k_ij;
    ptr[1] = m_l_ij;
    for (k=0; k<ptr[0]; k++) {
      for (l=0; l<ptr[1]; l++) {
        ptr[l*m_k_ij+k+2] = l*m_k_ij+k+idx;
      }
    }
    if (nv < nelems-1) {
      elem_buf = (void*)(((char*)elem_buf)+size);
    }
    nv++;
  }

  GP_Scatter(g_p, nv, subscripts, buf_ptr, buf_size, &size, 1);
  GP_Sync();

  /* Check contents of GP array */
  GP_Distribution(g_p, me, lo, hi);
  for (i=lo[0]; i<=hi[0]; i++) {
    ii = i - lo[0];
    subscript[0] = i;
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      subscript[1] = j;
      idx = j*N_I + i;
      GP_Access_element(g_p, subscript, &elem_buf, &size);
      ptr = (int*)elem_buf;
      m_k_ij = i%Q_I + 1;
      m_l_ij = j%Q_J + 1;
      if (size != sizeof(int)*(m_k_ij*m_l_ij+2)) {
        printf("p[%d] [%d,%d] Size actual: %d expected: %ld\n",
               me,i,j,size,(long)(sizeof(int)*(m_k_ij*m_l_ij+2)));
      }
      if (ptr[0] != m_k_ij) {
        printf("p[%d] [%d,%d] Dimension(3) i actual: %d expected: %d\n",me,i,j,ptr[0],m_k_ij);
      }
      if (ptr[1] != m_l_ij) {
        printf("p[%d] [%d,%d] Dimension(3) j actual: %d expected: %d\n",me,i,j,ptr[1],m_l_ij);
      }
      for (k=0; k<ptr[0]; k++) {
        for (l=0; l<ptr[1]; l++) {
          if (ptr[l*m_k_ij+k+2] != l*m_k_ij+k+idx) {
            printf("p[%d] Element(3) i: %d j: %d l: %d k: %d m_k_ij: %d idx: %d does not match: %d %d\n",
                me,i,j,l,k,m_k_ij,idx,ptr[l*ptr[0]+k+2],l*m_k_ij+k+idx);
          }
        }
      }
    }
  }
  if (me==0) printf("\nCompleted check of GP_Scatter\n");


  /* Clean up buffers and deallocate GP array */
  free(buf);
  free(buf_ptr);
  free(buf_size);
  free(subscripts);
  GP_Distribution(g_p, me, lo, hi);
  for (i=lo[0]; i<=hi[0]; i++) {
    subscript[0] = i;
    for (j=lo[1]; j<=hi[1]; j++) {
      subscript[1] = j;
      GP_Free(GP_Free_local_element(g_p, subscript));
    }
  }

  /* destroy Global Pointer array */
  GP_Destroy(g_p);

  /* Test set irregular distribution capability */
  dims3[0] = N_I;
  dims3[1] = N_J;
  dims3[2] = N_K;
  ndim = 3;

  g_p = GP_Create_handle();
  GP_Set_dimensions(g_p, ndim, dims3);

  /* get processor grid */
  idim = N_I;
  jdim = N_J;
  kdim = N_K;
  factor(nproc, idim, jdim, kdim, &pdi, &pdj, &pdk);

  /* construct map array */
  k = 0;
  mapc = (int*)malloc((pdi+pdj+pdk)*sizeof(int));
  for (i=0; i<pdi; i++) {
    mapc[k] = (i*idim)/pdi;
    k++;
  }
  for (i=0; i<pdj; i++) {
    mapc[k] = (i*jdim)/pdj;
    k++;
  }
  for (i=0; i<pdk; i++) {
    mapc[k] = (i*kdim)/pdk;
    k++;
  }
  blocks[0] = pdi;
  blocks[1] = pdj;
  blocks[2] = pdk;

  /* allocate array using an irregular distribution */
  GP_Set_irreg_distr(g_p, mapc, blocks);
  GP_Allocate(g_p);

  /* check local distribution */
  idx = me;
  k = idx%pdk;
  idx = (idx-k)/pdk;
  j = idx%pdj;
  i = (idx-j)/pdj;
  GP_Distribution(g_p, me, lo3, hi3);
  idx = 1;
  if (lo3[0] != (i*idim)/pdi) {
    printf("p[%d] Mismatch in irregular distribution array lo3[0]: %d expected: %d\n",
           me, lo3[0], (i*idim)/pdi);
    idx = 0;
  }
  if (hi3[0] != ((i+1)*idim)/pdi - 1) {
    printf("p[%d] Mismatch in irregular distribution array hi3[0]: %d expected: %d\n",
           me, hi3[0], ((i+1)*idim)/pdi-1);
    idx = 0;
  }
  if (lo3[1] != (j*jdim)/pdj) {
    printf("p[%d] Mismatch in irregular distribution array lo3[1]: %d expected: %d\n",
           me, lo3[1], (j*jdim)/pdj);
    idx = 0;
  }
  if (hi3[1] != ((j+1)*jdim)/pdj - 1) {
    printf("p[%d] Mismatch in irregular distribution array hi3[1]: %d expected: %d\n",
           me, hi3[1], ((j+1)*jdim)/pdj-1);
    idx = 0;
  }
  if (lo3[2] != (k*kdim)/pdk) {
    printf("p[%d] Mismatch in irregular distribution array lo3[2]: %d expected: %d\n",
           me, lo3[2], (k*kdim)/pdk);
    idx = 0;
  }
  if (hi3[2] != ((k+1)*kdim)/pdk - 1) {
    printf("p[%d] Mismatch in irregular distribution array hi3[2]: %d expected: %d\n",
           me, hi3[2], ((k+1)*kdim)/pdk-1);
    idx = 0;
  }

  GA_Igop(&idx,1,"*");
  if (idx == 1 && me == 0) {
    printf("\nCompleted check of GP_Set_irreg_distr\n");
  } else if (me == 0) {
    printf("\nFailed check of GP_Set_irreg_distr\n");
  }

  /* destroy Global Pointer array */
  GP_Destroy(g_p);

  /* Check set chunk capability */
  g_p = GP_Create_handle();
  GP_Set_dimensions(g_p, ndim, dims3);
  chunk[0] = -1;
  chunk[1] = -1;
  chunk[2] = dims3[2];
  GP_Set_chunk(g_p, chunk);
  GP_Allocate(g_p);

  GP_Distribution(g_p, me, lo3, hi3);
  idx = 1;
  if (lo3[2] != 0) {
    printf("p[%d] Mismatch in chunk distribution array lo3[2]: %d expected: %d\n",
           me, lo3[2], 0);
    idx = 0;
  }
  if (hi3[2] != dims3[2] - 1) {
    printf("p[%d] Mismatch in chunk distribution array hi3[2]: %d expected: %d\n",
           me, hi3[2], dims3[2]-1);
    idx = 0;
  }
  GA_Igop(&idx,1,"*");
  if (idx == 1 && me == 0) {
    printf("\nCompleted check of GP_Set_chunk\n");
  } else if (me == 0) {
    printf("\nFailed check of GP_Set_chunk\n");
  }

  /* destroy Global Pointer array */
  GP_Destroy(g_p);

}

int main(int argc, char **argv)
{
  int heap=20000, stack=20000;
  int me, nproc;

  MP_INIT(argc, argv);

  NGA_Initialize();
  GP_Initialize();
  me = NGA_Nodeid();
  nproc = NGA_Nnodes();

  if (me==0) {
    if (GA_Uses_fapi()) NGA_Error("Program runs with C array API only",0);
    printf("Using %ld processes\n", (long)nproc);
    fflush(stdout);
  }

  heap /= nproc;
  stack /= nproc;

  if (!MA_init(MT_F_DBL, stack, heap))
    NGA_Error("MA_init failed",stack+heap);

  do_work();

  if (me==0) printf("Terminating ..\n");
  GP_Terminate();
  NGA_Terminate();

  MP_FINALIZE();

  return 0;
}
