#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mock.h"
#include "abstract_ops.h"


void Mock_Abs_value(mock_ga_t *g_a)
{
  ITER_DECLARE_VARS(g_a)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)			\
      case GA_TYPE:					\
	{						\
	  ITER_DECLARE_BUFFER(g_a,C_TYPE)		\
	    ITER_INIT(g_a,C_TYPE)			\
	    ITER_BEGIN(g_a,C_TYPE)			\
	    assign_abs_##AT(*g_a_buf,*g_a_buf);		\
	  ITER_NEXT(g_a)				\
	    ITER_END					\
	    break;					\
	}
#include "types.xh"
#undef TYPE_CASE
    }
}


void Mock_Abs_value_patch(mock_ga_t *g_a, int *lo, int *hi)
{
  ITER_DECLARE_VARS_PATCH(g_a)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)			\
      case GA_TYPE:					\
	{						\
	  ITER_DECLARE_BUFFER(g_a,C_TYPE)		\
	    ITER_INIT_PATCH(g_a,C_TYPE,lo,hi)		\
	    ITER_BEGIN(g_a,C_TYPE)			\
	    assign_abs_##AT(*g_a_buf,*g_a_buf);		\
	  ITER_NEXT_PATCH(g_a)				\
	    ITER_END					\
	    break;					\
	}
#include "types.xh"
#undef TYPE_CASE
    }
}


void Mock_Access_block_grid(mock_ga_t *g_a, int index[], void *ptr, int ld[])
{

}


void Mock_Access_block(mock_ga_t *g_a, int idx, void *ptr, int ld[])
{

}


void Mock_Access_block_segment(mock_ga_t *g_a, int proc, void *ptr, int *len)
{

}


void Mock_Access_ghost_element(mock_ga_t *g_a, void *ptr, int subscript[], int ld[])
{

}


void Mock_Access_ghosts(mock_ga_t *g_a, int dims[], void *ptr, int ld[])
{

}


void Mock_Access(mock_ga_t *g_a, int lo[], int hi[], void *ptr, int ld[])
{

}


void Mock_Acc(mock_ga_t *g_a, int lo[], int hi[],void* buf, int ld[],void* alpha)
{

}


void Mock_Add_constant(mock_ga_t *g_a, void* alpha)
{
  ITER_DECLARE_VARS(g_a)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                \
      case GA_TYPE:				    \
	{					    \
	  C_TYPE value = *((C_TYPE*)alpha);	    \
	  ITER_DECLARE_BUFFER(g_a,C_TYPE)	    \
	    ITER_INIT(g_a,C_TYPE)		    \
	    ITER_BEGIN(g_a,C_TYPE)		    \
	    add_assign_##AT(*g_a_buf,value);	    \
	  ITER_NEXT(g_a)			    \
	    ITER_END				    \
	    break;				    \
	}
#include "types.xh"
#undef TYPE_CASE
    }
}

void Mock_Add_constant_patch(mock_ga_t *g_a, int *lo, int *hi,void *alpha)
{
  /*
  ITER_DECLARE_VARS_PATCH(g_a)
    
    switch (g_a->type){
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)			\
      case GA_TYPE:					\
	{						\
	  C_TYPE value = *((C_TYPE*)alpha);		\
	  ITER_DECLARE_BUFFER(g_a,C_TYPE)		\
	    ITER_INIT(g_a,C_TYPE)			\
	    ITER_BEGIN(g_a,C_TYPE)			\
	    add_assign_##AT(*g_a_buf,value);		\
	  ITIR_NEXT_PATCH(g_a)				\
	    ITER_END					\
	    break;					\
	}
#include "types.xh"
#undef TYPE_CASE
    }
  */
}

void Mock_Add_diagonal(mock_ga_t *g_a, int g_v)
{

}


void Mock_Add_patch(void * alpha, mock_ga_t *g_a, int alo[], int ahi[], void * beta, mock_ga_t *g_b, int blo[], int bhi[], mock_ga_t *g_c, int clo[], int chi[])
{

}


void Mock_Add(void *alpha, mock_ga_t *g_a, void* beta, mock_ga_t *g_b, mock_ga_t *g_c)
{

}


int Mock_Allocate(mock_ga_t *g_a)
{

}


int Mock_Assemble_duplicate(mock_ga_t *g_a, char *name, void *ptr)
{

}


void Mock_Brdcst(void *buf, int lenbuf, int root)
{

}


SingleComplex Mock_Cdot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


SingleComplex Mock_Cdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Cgemm(char ta, char tb, int m, int n, int k, SingleComplex alpha, mock_ga_t *g_a, mock_ga_t *g_b, SingleComplex beta, mock_ga_t *g_c )
{

}


void Mock_Cgop(SingleComplex x[], int n, char *op)
{

}


void Mock_Check_handle(mock_ga_t *g_a, char *string)
{

}


int Mock_Cluster_nnodes(void)
{

}


int Mock_Cluster_nodeid(void)
{

}


int Mock_Cluster_nprocs(int x)
{

}


int Mock_Cluster_procid(int x, int y)
{

}


int Mock_Cluster_proc_nodeid(int proc)
{

}


int Mock_Compare_distr(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


void Mock_Copy(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


void Mock_Copy_patch(char trans, mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[])
{

}


mock_ga_t* Mock_Create_config(int type, int ndim, int dims[], char *name, int chunk[], int p_handle)
{
    return Mock_Create(type, ndim, dims, name, chunk);
}


mock_ga_t* Mock_Create_ghosts_config(int type, int ndim, int dims[], int width[], char *name, int chunk[], int p_handle)
{
    return Mock_Create(type, ndim, dims, name, chunk);
}


mock_ga_t* Mock_Create_ghosts(int type, int ndim, int dims[], int width[], char *name, int chunk[])
{
    return Mock_Create(type, ndim, dims, name, chunk);
}


mock_ga_t* Mock_Create_ghosts_irreg_config(int type, int ndim, int dims[], int width[], char *name, int map[], int nblock[], int p_handle)
{
    return Mock_Create(type, ndim, dims, name, NULL);
}


mock_ga_t* Mock_Create_ghosts_irreg(int type, int ndim, int dims[], int width[], char *name, int map[], int nblock[])
{
    return Mock_Create(type, ndim, dims, name, NULL);
}


mock_ga_t* Mock_Create_handle(void)
{
    return (mock_ga_t*)malloc(sizeof(mock_ga_t));
}


mock_ga_t* Mock_Create(int type, int ndim, int dims[], char *name, int chunk[])
{
    int i;
    int nd_m1 = ndim-1;
    mock_ga_t *g_a = (mock_ga_t*)malloc(sizeof(mock_ga_t));

    g_a->type = type;
    g_a->ndim = ndim;
    g_a->nd_m1 = nd_m1;
    g_a->size = 1;
    g_a->itemsize = (int)type_to_size(type);
    memset(g_a->dims,        0, sizeof(int)*GA_MAX_DIM);
    memset(g_a->dims_m1,     0, sizeof(int)*GA_MAX_DIM);
    memset(g_a->strides,     0, sizeof(int)*GA_MAX_DIM);
    memset(g_a->backstrides, 0, sizeof(int)*GA_MAX_DIM);
    for (i=0; i<ndim; ++i) {
        g_a->dims[i] = dims[i];
        g_a->size *= dims[i];
    }
    /* setup important iterator variables */
    for (i=nd_m1; i>=0; --i) {
        g_a->dims_m1[i] = dims[i]-1;
        g_a->strides[i] = (i == nd_m1) ? 1 : g_a->strides[i+1]*g_a->dims[i+1];
        g_a->backstrides[i] = g_a->dims_m1[i]*g_a->strides[i];
    }
#if 0
    /* scale the back/strides based on itemsize */
    for (i=0; i<ndim; ++i) {
        g_a->strides[i] *= g_a->itemsize;
        g_a->backstrides[i] *= g_a->itemsize;
    }
#endif
    g_a->buf = malloc(g_a->size * g_a->itemsize);
#if 0
    printf("Mock_Create(%d, %d, dims, %s, chunk)\n", type, ndim, name);
    aprint("dims", g_a->dims, g_a->ndim);
    aprint("dims_m1", g_a->dims_m1, g_a->ndim);
    aprint("strides", g_a->strides, g_a->ndim);
    aprint("backstrides", g_a->backstrides, g_a->ndim);
#endif

    return g_a;
}


mock_ga_t* Mock_Create_irreg_config(int type, int ndim, int dims[],char *name, int block[], int map[], int p_handle)
{
    return Mock_Create(type, ndim, dims, name, NULL);
}


mock_ga_t* Mock_Create_irreg(int type, int ndim, int dims[],char *name, int block[], int map[])
{
    return Mock_Create(type, ndim, dims, name, NULL);
}


int Mock_Create_mutexes(int number)
{

}


double Mock_Ddot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


double Mock_Ddot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


int Mock_Deregister_type(int type)
{

}


void Mock_Destroy(mock_ga_t *g_a)
{

}


int Mock_Destroy_mutexes(void)
{

}


void Mock_Dgemm(char ta, char tb, int m, int n, int k, double alpha, mock_ga_t *g_a, mock_ga_t *g_b, double beta, mock_ga_t *g_c )
{

}


void Mock_Dgop(double x[], int n, char *op)
{

}


void Mock_Diag(mock_ga_t *g_a, int g_s, int g_v, void *eval)
{

}


void Mock_Diag_reuse(int reuse, mock_ga_t *g_a, int g_s, int g_v, void *eval)
{

}


void Mock_Diag_seq(mock_ga_t *g_a, int g_s, int g_v, void *eval)
{

}


void Mock_Diag_std(mock_ga_t *g_a, int g_v, void *eval)
{

}


void Mock_Diag_std_seq(mock_ga_t *g_a, int g_v, void *eval)
{

}


void Mock_Distribution(mock_ga_t *g_a, int iproc, int lo[], int hi[])
{

}


mock_ga_t* Mock_Duplicate(mock_ga_t *g_a, char* array_name)
{

}


void Mock_Elem_divide(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c)
{
    /* TODO FIXME: wasn't compiling */
#if 0
  ITER_DECLARE_VARS(g_a)
    ITER_DECLARE_VARS(g_b)
    ITER_DECLARE_VARS(g_c)
    
      switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                         \
        case GA_TYPE:                                        \
	  {							     \
	    ITER_DECLARE_BUFFER(g_a,C_TYPE)			     \
	      ITER_DECLARE_BUFFER(g_b,C_TYPE)			     \
	      ITER_DECLARE_BUFFER(g_c,C_TYPE)			     \
	      ITER_INIT(g_a,C_TYPE)				     \
	      ITER_INIT(g_b,C_TYPE)				     \
	      ITER_INIT(g_c,C_TYPE)				     \
	      ITER_BEGIN(g_a,C_TYPE)				     \
	      assign_div_##AT(*g_c_buf,*g_a_buf,*g_b_buf);	     \
	    ITER_NEXT(g_a)					     \
	      ITER_NEXT(g_b)					     \
	      ITER_NEXT(g_c)					     \
	      ITER_END						     \
	      break;						     \
	  }
#include "types.xh"
#undef TYPE_CASE
      }
#endif
}


void Mock_Elem_divide_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi,mock_ga_t *g_c, int *clo, int *chi)
{

}


void Mock_Elem_maximum(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c)
{
  ITER_DECLARE_VARS(g_a)
    ITER_DECLARE_VARS(g_b)
    ITER_DECLARE_VARS(g_c)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                         \
      case GA_TYPE:					     \
	{						     \
	    ITER_DECLARE_BUFFER(g_a,C_TYPE)		     \
	    ITER_DECLARE_BUFFER(g_b,C_TYPE)		     \
	    ITER_DECLARE_BUFFER(g_c,C_TYPE)		     \
	    ITER_INIT(g_a,C_TYPE)			     \
	    ITER_INIT(g_b,C_TYPE)			     \
	    ITER_INIT(g_c,C_TYPE)			     \
	    ITER_BEGIN(g_a,C_TYPE)			     \
	    assign_max_##AT(*g_c_buf,*g_a_buf,*g_b_buf);     \
	    ITER_NEXT(g_a)				     \
	    ITER_NEXT(g_b)				     \
	    ITER_NEXT(g_c)				     \
	    ITER_END					     \
	    break;					     \
	}
#include "types.xh"
#undef TYPE_CASE
    }

}


void Mock_Elem_maximum_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi,mock_ga_t *g_c, int *clo, int *chi)
{

}


void Mock_Elem_minimum(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c)
{
  ITER_DECLARE_VARS(g_a)
    ITER_DECLARE_VARS(g_b)
    ITER_DECLARE_VARS(g_c)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)				     \
      case GA_TYPE:						     \
        {							     \
          ITER_DECLARE_BUFFER(g_a,C_TYPE)                            \
	    ITER_DECLARE_BUFFER(g_b,C_TYPE)			     \
	    ITER_DECLARE_BUFFER(g_c,C_TYPE)			     \
	    ITER_INIT(g_a,C_TYPE)				     \
	    ITER_INIT(g_b,C_TYPE)				     \
	    ITER_INIT(g_c,C_TYPE)				     \
	    ITER_BEGIN(g_a,C_TYPE)				     \
            assign_max_##AT(*g_c_buf,*g_a_buf,*g_b_buf);	     \
          ITER_NEXT(g_a)					     \
	    ITER_NEXT(g_b)					     \
	    ITER_NEXT(g_c)					     \
	    ITER_END						     \
            break;						     \
        }
#include "types.xh"
#undef TYPE_CASE
    }

}


void Mock_Elem_minimum_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi,mock_ga_t *g_c, int *clo, int *chi)
{
  /*  ITER_DECLARE_VARS_PATCH(g_a)
  ITER_DECLARE_VARS_PATCH(g_b)
  ITER_DECLARE_VARS_PATCH(g_c)

    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                \
      case GA_TYPE:                               \
    {                                       \
      ITER_DECLARE_BUFFER(g_a,C_TYPE)             \
      ITER_DECLARE_BUFFER(g_b,C_TYPE)             \
      ITER_DECLARE_BUFFER(g_c,C_TYPE)             \
        ITER_BEGIN(g_a,C_TYPE,lo,hi)  \
        assign_abs_##AT(*g_a_buf,*g_a_buf); \
      ITER_NEXT_PATCH(g_a)                \
                ITER_END                            \
        break;                              \
    }
#include "types.xh"
#undef TYPE_CASE
    }
  */
}


void Mock_Elem_multiply(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c)
{
  ITER_DECLARE_VARS(g_a)
  ITER_DECLARE_VARS(g_b)
  ITER_DECLARE_VARS(g_c)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)               \
      case GA_TYPE:                                \
    {                                              \
      ITER_DECLARE_BUFFER(g_a,C_TYPE)              \
      ITER_DECLARE_BUFFER(g_b,C_TYPE)              \
      ITER_DECLARE_BUFFER(g_c,C_TYPE)              \
      ITER_INIT(g_a,C_TYPE)                        \
      ITER_INIT(g_b,C_TYPE)                        \
      ITER_INIT(g_c,C_TYPE)                        \
      ITER_BEGIN(g_a,C_TYPE)                       \
      assign_mul_##AT(*g_c_buf,*g_a_buf,*g_b_buf); \
      ITER_NEXT(g_a)                               \
      ITER_NEXT(g_b)                               \
      ITER_NEXT(g_c)                               \
      ITER_END                                     \
      break;                                       \
    }
#include "types.xh"
#undef TYPE_CASE
    }

}


void Mock_Elem_multiply_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi,mock_ga_t *g_c, int *clo, int *chi)
{

}


void Mock_Error(char *str, int code)
{

}


float Mock_Fdot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


float Mock_Fdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Fence(void)
{

}


void Mock_Fgop(float x[], int n, char *op)
{

}


void Mock_Fill(mock_ga_t *g_a, void *value)
{
  
  ITER_DECLARE_VARS(g_a)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                \
      case GA_TYPE:                                 \
        {                                           \
          ITER_DECLARE_BUFFER(g_a,C_TYPE)	    \
	    ITER_INIT(g_a,C_TYPE)		    \
	    ITER_BEGIN(g_a,C_TYPE)                  \
	    assign_##AT(*g_a_buf,*g_a_buf);	    \
          ITER_NEXT(g_a)                            \
            ITER_END                                \
            break;                                  \
        }
#include "types.xh"
#undef TYPE_CASE
    }      
}
  

void Mock_Fill_patch(mock_ga_t *g_a, int lo[], int hi[], void *val)
{

}


void Mock_Freemem(void* ptr)
{

}


void Mock_Gather_flat(mock_ga_t *g_a, void *v, int subsArray[], int n)
{

}


void Mock_Gather(mock_ga_t *g_a, void *v, int* subsArray[], int n)
{

}


void Mock_Get_block_info(mock_ga_t *g_a, int num_blocks[], int block_dims[])
{

}


int Mock_Get_debug(void)
{

}


void Mock_Get_diag(mock_ga_t *g_a, int g_v)
{

}


int Mock_Get_dimension(mock_ga_t *g_a)
{

}


void Mock_Get_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld)
{

}


void Mock_Get_ghost_block(mock_ga_t *g_a, int lo[], int hi[], void *buf, int ld[])
{

}


void Mock_Get(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[])
{

}


void* Mock_Getmem(int type, int nelem, int grp_id)
{

}


int Mock_Get_pgroup(mock_ga_t *g_a)
{

}


int Mock_Get_pgroup_size(int grp_id)
{

}


void Mock_Get_proc_grid(mock_ga_t *g_a, int dims[])
{

}


void Mock_Get_proc_index(mock_ga_t *g_a, int iproc, int subscript[])
{

}


void Mock_Gop(int type, void *x, int n, char *op)
{

}


int Mock_Has_ghosts(mock_ga_t *g_a)
{

}


int Mock_Idot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


int Mock_Idot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Igop(int x[], int n, char *op)
{

}


void Mock_Init_fence(void)
{

}


void Mock_Initialize_args(int *argc, char ***argv)
{

}


void Mock_Initialize_ltd(size_t limit)
{

}


void Mock_Initialize(void)
{

}


void Mock_Inquire(mock_ga_t *g_a, int *type, int *ndim, int dims[])
{

}


size_t Mock_Inquire_memory(void)
{

}


char* Mock_Inquire_name(mock_ga_t *g_a)
{

}


int Mock_Is_mirrored(mock_ga_t *g_a)
{

}


long Mock_Ldot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


long Mock_Ldot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Lgop(long x[], int n, char *op)
{

}


void Mock_List_nodeid(int *list, int nprocs)
{

}


long long Mock_Lldot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


long long Mock_Lldot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Llgop(long long x[], int n, char *op)
{

}


int Mock_Llt_solve(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


int Mock_Locate(mock_ga_t *g_a, int subscript[])
{

}


int Mock_Locate_nnodes(mock_ga_t *g_a, int lo[], int hi[])
{

}


int Mock_Locate_num_blocks(mock_ga_t *g_a, int lo[], int hi[])
{

}


int Mock_Locate_region(mock_ga_t *g_a, int lo[], int hi[], int map[], int procs[])
{

}


void Mock_Lock(int mutex)
{

}


void Mock_Lu_solve(char tran, mock_ga_t *g_a, mock_ga_t *g_b)
{

}


void Mock_Mask_sync(int first, int last)
{

}


void Mock_Matmul_patch_2d(char transa, char transb, void* alpha, void *beta, mock_ga_t *g_a, int ailo, int aihi, int ajlo, int ajhi, mock_ga_t *g_b, int bilo, int bihi, int bjlo, int bjhi, mock_ga_t *g_c, int cilo, int cihi, int cjlo, int cjhi)
{

}


void Mock_Matmul_patch(char transa, char transb, void* alpha, void *beta, mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[], mock_ga_t *g_c, int clo[], int chi[]) 
{

}


void Mock_Median(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c, int g_m)
{

}


void Mock_Median_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi, mock_ga_t *g_c, int *clo, int *chi, int g_m, int *mlo, int *mhi)
{

}


size_t Mock_Memory_avail(void)
{

}


int Mock_Memory_limited(void)
{

}


void Mock_Merge_distr_patch(mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[])
{

}


void Mock_Merge_mirrored(mock_ga_t *g_a)
{

}


void Mock_NbAcc(mock_ga_t *g_a, int lo[], int hi[],void* buf, int ld[],void* alpha, ga_nbhdl_t* nbhandle)
{

}


void Mock_Nbget_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize,void *buf, int *ld, ga_nbhdl_t *nbhandle)
{

}


void Mock_NbGet_ghost_dir(mock_ga_t *g_a, int mask[], ga_nbhdl_t* handle)
{

}


void Mock_NbGet(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle)
{

}


void Mock_Nblock(mock_ga_t *g_a, int *nblock)
{

}


void Mock_Nbput_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld, ga_nbhdl_t *nbhandle)
{

}


void Mock_NbPut(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle)
{

}


int Mock_NbTest(ga_nbhdl_t* nbhandle)
{

}


void Mock_NbWait(ga_nbhdl_t* nbhandle)
{

}


int Mock_Ndim(mock_ga_t *g_a)
{

}


int Mock_Nnodes(void)
{

}


int Mock_Nodeid(void)
{

}


void Mock_Norm1(mock_ga_t *g_a, double *nm)
{

}


void Mock_Norm_infinity(mock_ga_t *g_a, double *nm)
{

}


void Mock_Periodic_acc(mock_ga_t *g_a, int lo[], int hi[],void* buf, int ld[],void* alpha)
{

}


void Mock_Periodic_get(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[])
{

}


void Mock_Periodic_put(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[])
{

}


int Mock_Pgroup_absolute_id(int pgroup, int pid)
{

}


void Mock_Pgroup_brdcst(int grp, void *buf, int lenbuf, int root)
{

}


void Mock_Pgroup_cgop(int grp, SingleComplex x[], int n, char *op)
{

}


int Mock_Pgroup_create(int *list, int count)
{

}


int Mock_Pgroup_destroy(int grp)
{

}


void Mock_Pgroup_dgop(int grp, double x[], int n, char *op)
{

}


void Mock_Pgroup_fgop(int grp, float x[], int n, char *op)
{

}


int Mock_Pgroup_get_default(void)
{

}


int Mock_Pgroup_get_mirror(void)
{

}


int Mock_Pgroup_get_world(void)
{

}


void Mock_Pgroup_igop(int grp, int x[], int n, char *op)
{

}


void Mock_Pgroup_lgop(int grp, long x[], int n, char *op)
{

}


void Mock_Pgroup_llgop(int grp, long long x[], int n, char *op)
{

}


int Mock_Pgroup_nnodes(int grp_id)
{

}


int Mock_Pgroup_nodeid(int grp_id)
{

}


void Mock_Pgroup_set_default(int p_handle)
{

}


int Mock_Pgroup_split(int grp_id, int num_group)
{

}


int Mock_Pgroup_split_irreg(int grp_id, int color)
{

}


void Mock_Pgroup_sync(int grp_id)
{

}


void Mock_Pgroup_zgop(int grp, DoubleComplex x[], int n, char *op)
{

}


void Mock_Print_distribution(mock_ga_t *g_a)
{

}


void Mock_Print_file(FILE *file, mock_ga_t *g_a)
{

}


/* just for debugging purposes */
void Mock_Print(mock_ga_t *g_a)
{
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                        \
        case GA_TYPE:                                       \
            {                                               \
                int i,j;                                    \
                int coords[GA_MAX_DIM] = {0};               \
                C_TYPE *ptr = (C_TYPE*)(g_a->buf);          \
                for (i=0; i<g_a->size; ++i) {               \
                    printf("%d (%d", i, coords[0]);         \
                    for (j=1; j<g_a->ndim; ++j) {           \
                        printf(",%d", coords[j]);           \
                    }                                       \
                    printf(") " AT "\n", ptr[i]);           \
                    for (j=g_a->nd_m1; j>=0; --j) {         \
                        if (coords[j] < g_a->dims_m1[j]) {  \
                            ++coords[j];                    \
                            break;                          \
                        }                                   \
                        else {                              \
                            coords[j] = 0;                  \
                        }                                   \
                    }                                       \
                }                                           \
                break;                                      \
            }
        TYPE_CASE(C_INT,int,"%d")
        TYPE_CASE(C_LONG,long,"%ld")
        TYPE_CASE(C_LONGLONG,long long,"%lld")
        TYPE_CASE(C_FLOAT,float,"%f")
        TYPE_CASE(C_DBL,double,"%f")
#undef TYPE_CASE
#define TYPE_CASE(GA_TYPE,C_TYPE)                                   \
        case GA_TYPE:                                               \
            {                                                       \
                int i,j;                                            \
                int coords[GA_MAX_DIM] = {0};                       \
                C_TYPE *ptr = (C_TYPE*)(g_a->buf);                  \
                for (i=0; i<g_a->size; ++i) {                       \
                    printf("%d (%d", i, coords[0]);                 \
                    for (j=1; j<g_a->ndim; ++j) {                   \
                        printf(",%d", coords[j]);                   \
                    }                                               \
                    printf(") %f+%fi\n", ptr[i].real, ptr[i].imag); \
                    for (j=g_a->nd_m1; j>=0; --j) {                 \
                        if (coords[j] <= g_a->dims_m1[j]) {         \
                            ++coords[j];                            \
                            break;                                  \
                        }                                           \
                        else {                                      \
                            coords[j] = 0;                          \
                        }                                           \
                    }                                               \
                }                                                   \
                break;                                              \
            }
        TYPE_CASE(C_SCPL,SingleComplex)
        TYPE_CASE(C_DCPL,DoubleComplex)
#undef TYPE_CASE
    }
}


void Mock_Print_patch_2d(mock_ga_t *g_a, int ilo, int ihi, int jlo, int jhi, int pretty)
{

}


void Mock_Print_patch(mock_ga_t *g_a, int lo[], int hi[], int pretty)
{

}


void Mock_Print_stats(void)
{

}


void Mock_Proc_topology(mock_ga_t *g_a, int proc, int coord[])
{

}


void Mock_Put_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld)
{

}


void Mock_Put(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[])
{

}


void Mock_Randomize(mock_ga_t *g_a, void *value)
{
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                                 \
        case GA_TYPE:                                           \
            {                                              \
                int i;                                     \
                C_TYPE *ptr = (C_TYPE*)(g_a->buf);                   \
                C_TYPE val = *((C_TYPE*)value);                      \
                for (i=0; i<g_a->size; ++i) {              \
                    assign_rand_##AT(ptr[i],val);          \
                }                                          \
                break;                                     \
            }
#include "types.xh"
#undef TYPE_CASE
    }

}


long Mock_Read_inc(mock_ga_t *g_a, int subscript[], long inc)
{

}


void Mock_Recip(mock_ga_t *g_a)
{

}


void Mock_Recip_patch(mock_ga_t *g_a, int *lo, int *hi)
{

}


void Mock_Register_stack_memory(void * (*ext_alloc)(size_t, int, char *), void (*ext_free)(void *))
{

}


int Mock_Register_type(size_t bytes)
{

}


void Mock_Release_block_grid(mock_ga_t *g_a, int index[])
{

}


void Mock_Release_block(mock_ga_t *g_a, int idx)
{

}


void Mock_Release_block_segment(mock_ga_t *g_a, int idx)
{

}


void Mock_Release_ghost_element(mock_ga_t *g_a, int index[])
{

}


void Mock_Release_ghosts(mock_ga_t *g_a)
{

}


void Mock_Release(mock_ga_t *g_a, int lo[], int hi[])
{

}


void Mock_Release_update_block_grid(mock_ga_t *g_a, int index[])
{

}


void Mock_Release_update_block(mock_ga_t *g_a, int idx)
{

}


void Mock_Release_update_block_segment(mock_ga_t *g_a, int idx)
{

}


void Mock_Release_update_ghost_element(mock_ga_t *g_a, int index[])
{

}


void Mock_Release_update_ghosts(mock_ga_t *g_a)
{

}


void Mock_Release_update(mock_ga_t *g_a, int lo[], int hi[])
{

}


void Mock_Scale_cols(mock_ga_t *g_a, int g_v)
{

}


void Mock_Scale(mock_ga_t *g_a, void *value)
{
  /*
    ITER_DECLARE_VARS(g_a)
    
    switch (g_a->type) {
    #define TYPE_CASE(GA_TYPE,C_TYPE,AT)	       \
    case GA_TYPE:				       \
    {						       \
    C_TYPE value = *((C_TYPE*)alpha);		       \
    ITER_DECLARE_BUFFER(g_a,C_TYPE)		       \
    ITER_INIT(g_a,C_TYPE)			       \
    ITER_BEGIN(g_a,C_TYPE)			       \
    add_assign_##AT(*g_a_buf,value);		       \
    ITER_NEXT(g_a)				       \
    ITER_END					       \
    break;					       \
        }
	#include "types.xh"
	#undef TYPE_CASE
	}
  */
}


void Mock_Scale_patch(mock_ga_t *g_a, int lo[], int hi[], void *alpha)
{

}


void Mock_Scale_rows(mock_ga_t *g_a, int g_v)
{

}


void Mock_Scan_add(mock_ga_t *g_a, mock_ga_t *g_b, int g_sbit, int lo, int hi, int excl)
{

}


void Mock_Scan_copy(mock_ga_t *g_a, mock_ga_t *g_b, int g_sbit, int lo, int hi)
{

}


void Mock_Scatter_acc_flat(mock_ga_t *g_a, void *v, int subsArray[], int n, void *alpha)
{

}


void Mock_Scatter_acc(mock_ga_t *g_a, void *v, int* subsArray[], int n, void *alpha)
{

}


void Mock_Scatter_flat(mock_ga_t *g_a, void *v, int subsArray[], int n)
{

}


void Mock_Scatter(mock_ga_t *g_a, void *v, int* subsArray[], int n)
{

}


void Mock_Select_elem(mock_ga_t *g_a, char* op, void* val, int *index)
{

}


void Mock_Set_array_name(mock_ga_t *g_a, char *name)
{

}


void Mock_Set_block_cyclic(mock_ga_t *g_a, int dims[])
{

}


void Mock_Set_block_cyclic_proc_grid(mock_ga_t *g_a, int block[], int proc_grid[])
{

}


void Mock_Set_chunk(mock_ga_t *g_a, int chunk[])
{

}


void Mock_Set_data(mock_ga_t *g_a, int ndim, int dims[], int type)
{

}


void Mock_Set_debug(int flag)
{

}


void Mock_Set_diagonal(mock_ga_t *g_a, int g_v)
{

}


void Mock_Set_ghost_corner_flag(mock_ga_t *g_a, int flag)
{

}


void Mock_Set_ghosts(mock_ga_t *g_a, int width[])
{

}


void Mock_Set_irreg_distr(mock_ga_t *g_a, int map[], int block[])
{

}


void Mock_Set_irreg_flag(mock_ga_t *g_a, int flag)
{

}


void Mock_Set_memory_limit(size_t limit)
{

}


void Mock_Set_pgroup(mock_ga_t *g_a, int p_handle)
{

}


void Mock_Set_restricted(mock_ga_t *g_a, int list[], int size)
{

}


void Mock_Set_restricted_range(mock_ga_t *g_a, int lo_proc, int hi_proc)
{

}


void Mock_Sgemm(char ta, char tb, int m, int n, int k, float alpha, mock_ga_t *g_a, mock_ga_t *g_b, float beta, mock_ga_t *g_c )
{

}


void Mock_Shift_diagonal(mock_ga_t *g_a, void *c)
{

}


int Mock_Solve(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


int Mock_Spd_invert(mock_ga_t *g_a)
{

}


void Mock_Step_bound_info(int g_xx, int g_vv, int g_xxll, int g_xxuu, void *boundmin, void *wolfemin, void *boundmax)
{

}


void Mock_Step_bound_info_patch(int g_xx, int *xxlo, int *xxhi, int g_vv, int *vvlo, int *vvhi, int g_xxll, int *xxlllo, int *xxllhi, int g_xxuu, int *xxuulo, int *xxuuhi, void *boundmin, void *wolfemin, void *boundmax)
{

}


void Mock_Step_max(mock_ga_t *g_a, mock_ga_t *g_b, void *step)
{

}


void Mock_Step_max_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi, void *step)
{

}


void Mock_Strided_acc(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[], void *alpha)
{

}


void Mock_Strided_get(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[])
{

}


void Mock_Strided_put(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[])
{

}


void Mock_Summarize(int verbose)
{

}


void Mock_Symmetrize(mock_ga_t *g_a)
{

}


void Mock_Sync(void)
{

}


void Mock_Terminate(void)
{

}


int Mock_Total_blocks(mock_ga_t *g_a)
{

}


void Mock_Transpose(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


void Mock_Unlock(int mutex)
{

}


int Mock_Update_ghost_dir(mock_ga_t *g_a, int dimension, int idir, int flag)
{

}


void Mock_Update_ghosts(mock_ga_t *g_a)
{

}


int Mock_Uses_fapi(void)
{

}


int Mock_Uses_ma(void)
{

}


int Mock_Uses_proc_grid(mock_ga_t *g_a)
{

}


int Mock_Valid_handle(mock_ga_t *g_a)
{
    return g_a != NULL;
}


int Mock_Verify_handle(mock_ga_t *g_a)
{

}


double Mock_Wtime(void)
{

}


DoubleComplex Mock_Zdot(mock_ga_t *g_a, mock_ga_t *g_b)
{

}


DoubleComplex Mock_Zdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[])
{

}


void Mock_Zero_diagonal(mock_ga_t *g_a)
{

}


void Mock_Zero(mock_ga_t *g_a)
{
  
  ITER_DECLARE_VARS(g_a)
    
    switch (g_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)			    \
      case GA_TYPE:					    \
        {						    \
          ITER_DECLARE_BUFFER(g_a,C_TYPE)                   \
	    ITER_INIT(g_a,C_TYPE)			    \
	    ITER_BEGIN(g_a,C_TYPE)			    \
            assign_zero_##AT(*g_a_buf);			    \
          ITER_NEXT(g_a)				    \
            ITER_END					    \
            break;					    \
        }
#include "types.xh"
#undef TYPE_CASE
    }     
}


void Mock_Zero_patch(mock_ga_t *g_a, int lo[], int hi[])
{

}


void Mock_Zgemm(char ta, char tb, int m, int n, int k, DoubleComplex alpha, mock_ga_t *g_a, mock_ga_t *g_b, DoubleComplex beta, mock_ga_t *g_c )
{

}


void Mock_Zgop(DoubleComplex x[], int n, char *op)
{

}
