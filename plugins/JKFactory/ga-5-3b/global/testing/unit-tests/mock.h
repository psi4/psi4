#ifndef GLOBAL_TESTING_UNIT_TEST_MOCK_H
#define GLOBAL_TESTING_UNIT_TEST_MOCK_H

#include <assert.h>
#include <stdio.h>

#include "ga.h"
#include "typesf2c.h"
#include "abstract_ops.h"

typedef struct Mock_GA
{
    int type;
    int ndim;
    int dims[GA_MAX_DIM];
    int dims_m1[GA_MAX_DIM]; /* for convenience */
    int strides[GA_MAX_DIM]; /* for convenience */
    int backstrides[GA_MAX_DIM]; /* for convenience */
    int size; /* for convenience, total number of elements */
    int itemsize; /* for convenience, sizeof(type) */
    int nd_m1; /* for convenience, ndim-1 */
    void *buf;
} mock_ga_t;

static size_t type_to_size(int type)
{
    size_t size;
    switch (type) {
        case C_INT:      size = sizeof(int); break;
        case C_LONG:     size = sizeof(long); break;
        case C_LONGLONG: size = sizeof(long long); break;
        case C_FLOAT:    size = sizeof(float); break;
        case C_DBL:      size = sizeof(double); break;
        case C_SCPL:     size = sizeof(SingleComplex); break;
        case C_DCPL:     size = sizeof(DoubleComplex); break;
        default: assert(0);
    }
    return size;
}

/* all the variables we need for the iterator, prefixed with the given name */
#define ITER_DECLARE_VARS(the_ga)       \
char *the_ga##_bytebuf = (the_ga)->buf; \
int the_ga##_size = (the_ga)->size;     \
int the_ga##_i;

/* all the variables we need for the iterator, prefixed with the given name
 * Note: use this when iterating over a lo/hi patch */
#define ITER_DECLARE_VARS_PATCH(the_ga) \
char *the_ga##_bytebuf = (the_ga)->buf; \
int the_ga##_ndim = (the_ga)->ndim;     \
int the_ga##_nd_m1 = (the_ga)->nd_m1;   \
int the_ga##_coords[GA_MAX_DIM];        \
int the_ga##_lo[GA_MAX_DIM];            \
int the_ga##_hi[GA_MAX_DIM];            \
int the_ga##_size = 1;                  \
int the_ga##_i;                         \
int the_ga##_j;


/* to define the user-visible buffer on which to perform the user-defined op */
#define ITER_DECLARE_BUFFER(the_ga,TYPE) TYPE *the_ga##_buf = NULL;

/* setup the user-visible buffer and the iterator variables 
 * Note: use this when iterating over the entire array buffer i.e. no lo/hi */
#define ITER_CAST_BUFFER(the_ga,TYPE) the_ga##_buf = (TYPE*)the_ga##_bytebuf;

/* setup the user-visible buffer and the iterator variables 
 * Note: use this when iterating over the entire array buffer i.e. no lo/hi */
#define ITER_BEGIN(the_ga,TYPE) \
for (the_ga##_i = 0; the_ga##_i<the_ga##_size; ++the_ga##_i) {

#define ITER_INIT(the_ga,TYPE) \
ITER_CAST_BUFFER(the_ga,TYPE)

/* setup the user-visible buffer and the iterator variables 
 * Note: use this when iterating over a lo/hi patch */
#define ITER_INIT_PATCH(the_ga,TYPE,lo,hi)                              \
for (the_ga##_i = 0; the_ga##_i < the_ga##_ndim; ++ the_ga##_i) {       \
    the_ga##_coords[the_ga##_i] = 0;                                    \
    the_ga##_lo[the_ga##_i] = (lo)[the_ga##_i];                         \
    the_ga##_hi[the_ga##_i] = (hi)[the_ga##_i];                         \
    the_ga##_size *= the_ga##_hi[the_ga##_i]-the_ga##_lo[the_ga##_i]+1; \
    the_ga##_bytebuf += (the_ga)->strides[the_ga##_i] * the_ga##_lo[the_ga##_i] * (the_ga)->itemsize;                                                   \
}                                                                       \
ITER_CAST_BUFFER(the_ga,TYPE)

/* iterate over the buffer */
#define ITER_NEXT(the_ga) ++the_ga##_buf;

/* iterate over the buffer */
#define ITER_NEXT_PATCH(the_ga)                                       \
    for (the_ga##_j=the_ga##_nd_m1; the_ga##_j>=0; --the_ga##_j) {    \
        if (the_ga##_coords[the_ga##_j] <= the_ga##_hi[the_ga##_j]) { \
            ++the_ga##_coords[the_ga##_j];                            \
            the_ga##_buf += (the_ga)->strides[the_ga##_j];            \
            break;                                                    \
        }                                                             \
        else {                                                        \
            the_ga##_coords[the_ga##_j] = the_ga##_lo[the_ga##_j];    \
            the_ga##_buf -= (the_ga)->backstrides[the_ga##_j];        \
        }                                                             \
    }

/* end of loop */
#define ITER_END }

static void mock_data(mock_ga_t *mock_a, int g_a)
{
    ITER_DECLARE_VARS(mock_a)

    switch(mock_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                    \
        case GA_TYPE:                                   \
            {                                           \
                ITER_DECLARE_BUFFER(mock_a,C_TYPE)      \
                ITER_INIT(mock_a,C_TYPE)                \
                ITER_BEGIN(mock_a,C_TYPE)               \
                if (mock_a_i%2) {                       \
                    assign_##AT(*mock_a_buf, mock_a_i); \
                }                                       \
                else {                                  \
                    assign_##AT(*mock_a_buf,-mock_a_i); \
                }                                       \
                ITER_NEXT(mock_a)                       \
                ITER_END                                \
                break;                                  \
            }
        TYPE_CASE(C_INT,int,reg)
        TYPE_CASE(C_LONG,long,reg)
        TYPE_CASE(C_LONGLONG,long long,reg)
        TYPE_CASE(C_FLOAT,float,reg)
        TYPE_CASE(C_DBL,double,reg)
#undef TYPE_CASE
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                    \
        case GA_TYPE:                                   \
            {                                           \
                ITER_DECLARE_BUFFER(mock_a,C_TYPE)      \
                ITER_INIT(mock_a,C_TYPE)                \
                ITER_BEGIN(mock_a,C_TYPE)               \
                C_TYPE value;                           \
                if (mock_a_i%2) {                       \
                    value.real =  mock_a_i;             \
                    value.imag = -mock_a_i;             \
                }                                       \
                else {                                  \
                    value.real = -mock_a_i;             \
                    value.imag =  mock_a_i;             \
                }                                       \
                assign_##AT(*mock_a_buf,value);         \
                ITER_NEXT(mock_a)                       \
                ITER_END                                \
                break;                                  \
            }
        TYPE_CASE(C_SCPL,SingleComplex,cpl)
        TYPE_CASE(C_DCPL,DoubleComplex,cpl)
#undef TYPE_CASE
    }
}

static void mock_to_global(mock_ga_t *mock_a, int g_a)
{
    if (0 == GA_Nodeid()) {
        int i;
        int lo[GA_MAX_DIM];
        int hi[GA_MAX_DIM];
        int ld[GA_MAX_DIM];
        for (i=0; i<mock_a->ndim; ++i) {
            lo[i] = 0;
            hi[i] = mock_a->dims[i]-1;
        }
        for (i=0; i<mock_a->ndim-1; ++i) {
            ld[i] = hi[i+1]-lo[i+1]+1;
        }
        NGA_Put(g_a, lo, hi, mock_a->buf, ld);
    }
    GA_Sync();
}

static void global_to_mock(int g_a, mock_ga_t *mock_a)
{
    int i;
    int lo[GA_MAX_DIM];
    int hi[GA_MAX_DIM];
    int ld[GA_MAX_DIM];
    for (i=0; i<mock_a->ndim; ++i) {
        lo[i] = 0;
        hi[i] = mock_a->dims[i]-1;
    }
    for (i=0; i<mock_a->ndim-1; ++i) {
        ld[i] = hi[i+1]-lo[i+1]+1;
    }
    GA_Sync();
    NGA_Get(g_a, lo, hi, mock_a->buf, ld);
    GA_Sync();
}

static void print_val(void *val, int type)
{
    switch(type) {
        case C_INT:      printf("%d",     *((int*)val)); break;
        case C_LONG:     printf("%ld",    *((long*)val)); break;
        case C_LONGLONG: printf("%lld",   *((long long*)val)); break;
        case C_FLOAT:    printf("%f",     *((float*)val)); break;
        case C_DBL:      printf("%f",     *((double*)val)); break;
        case C_SCPL:     printf("%f+%fi", ((SingleComplex*)val)->real,((SingleComplex*)val)->imag); break;
        case C_DCPL:     printf("%f+%fi", ((DoubleComplex*)val)->real,((DoubleComplex*)val)->imag); break;
        default: assert(0);
    }
}

static int neq_mock(mock_ga_t *mock_a, mock_ga_t *mock_b, int *idx)
{
    ITER_DECLARE_VARS(mock_a)
    ITER_DECLARE_VARS(mock_b)

    *idx=-1;
    switch(mock_a->type) {
#define TYPE_CASE(GA_TYPE,C_TYPE,AT)                      \
        case GA_TYPE:                                     \
            {                                             \
                ITER_DECLARE_BUFFER(mock_a,C_TYPE)        \
                ITER_DECLARE_BUFFER(mock_b,C_TYPE)        \
                ITER_INIT(mock_a,C_TYPE)                  \
                ITER_INIT(mock_b,C_TYPE)                  \
                ITER_BEGIN(mock_a,C_TYPE)                 \
                if (!eq_##AT(*mock_a_buf, *mock_b_buf)) { \
                    *idx = mock_a_i;                      \
                    return 1;                             \
                }                                         \
                ITER_NEXT(mock_a)                         \
                ITER_NEXT(mock_b)                         \
                ITER_END                                  \
                break;                                    \
            }
#include "types.xh"
#undef TYPE_CASE
    }
    return 0;
}

void          Mock_Abs_value(mock_ga_t *g_a); 
void          Mock_Abs_value_patch(mock_ga_t *g_a, int *lo, int *hi);
void          Mock_Access_block_grid(mock_ga_t *g_a, int index[], void *ptr, int ld[]);
void          Mock_Access_block(mock_ga_t *g_a, int idx, void *ptr, int ld[]);
void          Mock_Access_block_segment(mock_ga_t *g_a, int proc, void *ptr, int *len);
void          Mock_Access_ghost_element(mock_ga_t *g_a,  void *ptr, int subscript[], int ld[]);
void          Mock_Access_ghosts(mock_ga_t *g_a, int dims[], void *ptr, int ld[]);
void          Mock_Access(mock_ga_t *g_a, int lo[], int hi[], void *ptr, int ld[]);
void          Mock_Acc(mock_ga_t *g_a, int lo[], int hi[],void* buf,int ld[],void* alpha);
void          Mock_Add_constant(mock_ga_t *g_a, void* alpha);
void          Mock_Add_constant_patch(mock_ga_t *g_a,int *lo,int *hi,void *alpha);
void          Mock_Add_diagonal(mock_ga_t *g_a, int g_v);
void          Mock_Add_patch(void * alpha, mock_ga_t *g_a, int alo[], int ahi[], void * beta,  mock_ga_t *g_b, int blo[], int bhi[], mock_ga_t *g_c, int clo[], int chi[]);
void          Mock_Add(void *alpha, mock_ga_t *g_a, void* beta, mock_ga_t *g_b, mock_ga_t *g_c); 
int           Mock_Allocate(mock_ga_t *g_a);
int           Mock_Assemble_duplicate(mock_ga_t *g_a, char *name, void *ptr);
void          Mock_Brdcst(void *buf, int lenbuf, int root);
SingleComplex Mock_Cdot(mock_ga_t *g_a, mock_ga_t *g_b); 
SingleComplex Mock_Cdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Cgemm(char ta, char tb, int m, int n, int k, SingleComplex alpha, mock_ga_t *g_a, mock_ga_t *g_b, SingleComplex beta, mock_ga_t *g_c );
void          Mock_Cgop(SingleComplex x[], int n, char *op);
void          Mock_Check_handle(mock_ga_t *g_a, char *string);
int           Mock_Cluster_nnodes(void);
int           Mock_Cluster_nodeid(void);
int           Mock_Cluster_nprocs(int x);
int           Mock_Cluster_procid(int x, int y);
int           Mock_Cluster_proc_nodeid(int proc);
int           Mock_Compare_distr(mock_ga_t *g_a, mock_ga_t *g_b); 
void          Mock_Copy(mock_ga_t *g_a, mock_ga_t *g_b); 
void          Mock_Copy_patch(char trans, mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[]);
mock_ga_t *   Mock_Create_config(int type,int ndim,int dims[], char *name, int chunk[], int p_handle);
mock_ga_t *   Mock_Create_ghosts_config(int type,int ndim,int dims[], int width[], char *name, int chunk[], int p_handle);
mock_ga_t *   Mock_Create_ghosts(int type,int ndim,int dims[], int width[], char *name, int chunk[]);
mock_ga_t *   Mock_Create_ghosts_irreg_config(int type,int ndim,int dims[], int width[], char *name, int map[], int nblock[], int p_handle);
mock_ga_t *   Mock_Create_ghosts_irreg(int type,int ndim,int dims[], int width[], char *name, int map[], int nblock[]);
mock_ga_t *   Mock_Create_handle(void);
mock_ga_t *   Mock_Create(int type,int ndim,int dims[], char *name, int chunk[]);
mock_ga_t *   Mock_Create_irreg_config(int type,int ndim,int dims[],char *name, int block[], int map[], int p_handle);
mock_ga_t *   Mock_Create_irreg(int type,int ndim,int dims[],char *name, int block[], int map[]);
int           Mock_Create_mutexes(int number);
double        Mock_Ddot(mock_ga_t *g_a, mock_ga_t *g_b); 
double        Mock_Ddot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
int           Mock_Deregister_type(int type);
void          Mock_Destroy(mock_ga_t *g_a);
int           Mock_Destroy_mutexes(void);
void          Mock_Dgemm(char ta, char tb, int m, int n, int k, double alpha, mock_ga_t *g_a, mock_ga_t *g_b, double beta, mock_ga_t *g_c );
void          Mock_Dgop(double x[], int n, char *op);
void          Mock_Diag(mock_ga_t *g_a, int g_s, int g_v, void *eval);
void          Mock_Diag_reuse(int reuse, mock_ga_t *g_a, int g_s, int g_v, void *eval);
void          Mock_Diag_seq(mock_ga_t *g_a, int g_s, int g_v, void *eval);
void          Mock_Diag_std(mock_ga_t *g_a, int g_v, void *eval);
void          Mock_Diag_std_seq(mock_ga_t *g_a, int g_v, void *eval);
void          Mock_Distribution(mock_ga_t *g_a, int iproc, int lo[], int hi[]); 
mock_ga_t *   Mock_Duplicate(mock_ga_t *g_a, char* array_name);
void          Mock_Elem_divide(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c);
void          Mock_Elem_divide_patch(mock_ga_t *g_a,int *alo,int *ahi, mock_ga_t *g_b,int *blo,int *bhi,mock_ga_t *g_c,int *clo,int *chi);
void          Mock_Elem_maximum(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c);
void          Mock_Elem_maximum_patch(mock_ga_t *g_a,int *alo,int *ahi, mock_ga_t *g_b,int *blo,int *bhi,mock_ga_t *g_c,int *clo,int *chi);
void          Mock_Elem_minimum(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c);
void          Mock_Elem_minimum_patch(mock_ga_t *g_a,int *alo,int *ahi, mock_ga_t *g_b,int *blo,int *bhi,mock_ga_t *g_c,int *clo,int *chi);
void          Mock_Elem_multiply(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c);
void          Mock_Elem_multiply_patch(mock_ga_t *g_a,int *alo,int *ahi, mock_ga_t *g_b,int *blo,int *bhi,mock_ga_t *g_c,int *clo,int *chi);
void          Mock_Error(char *str, int code);
float         Mock_Fdot(mock_ga_t *g_a, mock_ga_t *g_b);
float         Mock_Fdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Fence(void);
void          Mock_Fgop(float x[], int n, char *op);
void          Mock_Fill(mock_ga_t *g_a, void *value);
void          Mock_Fill_patch(mock_ga_t *g_a, int lo[], int hi[], void *val);
void          Mock_Freemem(void* ptr);
void          Mock_Gather_flat(mock_ga_t *g_a, void *v, int subsArray[], int n);
void          Mock_Gather(mock_ga_t *g_a, void *v, int* subsArray[], int n);
void          Mock_Get_block_info(mock_ga_t *g_a, int num_blocks[], int block_dims[]);
int           Mock_Get_debug(void);
void          Mock_Get_diag(mock_ga_t *g_a, int g_v);
int           Mock_Get_dimension(mock_ga_t *g_a);
void          Mock_Get_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld);
void          Mock_Get_ghost_block(mock_ga_t *g_a, int lo[], int hi[], void *buf, int ld[]);
void          Mock_Get(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[]); 
void*         Mock_Getmem(int type, int nelem, int grp_id);
int           Mock_Get_pgroup(mock_ga_t *g_a);
int           Mock_Get_pgroup_size(int grp_id);
void          Mock_Get_proc_grid(mock_ga_t *g_a, int dims[]);
void          Mock_Get_proc_index(mock_ga_t *g_a, int iproc, int subscript[]);
void          Mock_Gop(int type, void *x, int n, char *op);
int           Mock_Has_ghosts(mock_ga_t *g_a);
int           Mock_Idot(mock_ga_t *g_a, mock_ga_t *g_b);
int           Mock_Idot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Igop(int x[], int n, char *op);
void          Mock_Init_fence(void);
void          Mock_Initialize_args(int *argc, char ***argv);
void          Mock_Initialize_ltd(size_t limit);
void          Mock_Initialize(void);
void          Mock_Inquire(mock_ga_t *g_a, int *type, int *ndim, int dims[]);
size_t        Mock_Inquire_memory(void);
char*         Mock_Inquire_name(mock_ga_t *g_a);
int           Mock_Is_mirrored(mock_ga_t *g_a);
long          Mock_Ldot(mock_ga_t *g_a, mock_ga_t *g_b);
long          Mock_Ldot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Lgop(long x[], int n, char *op);
void          Mock_List_nodeid(int *list, int nprocs);
long long     Mock_Lldot(mock_ga_t *g_a, mock_ga_t *g_b);
long long     Mock_Lldot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Llgop(long long x[], int n, char *op);
int           Mock_Llt_solve(mock_ga_t *g_a, mock_ga_t *g_b);
int           Mock_Locate(mock_ga_t *g_a, int subscript[]);
int           Mock_Locate_nnodes(mock_ga_t *g_a, int lo[], int hi[]);
int           Mock_Locate_num_blocks(mock_ga_t *g_a, int lo[], int hi[]);
int           Mock_Locate_region(mock_ga_t *g_a,int lo[],int hi[],int map[],int procs[]);
void          Mock_Lock(int mutex);
void          Mock_Lu_solve(char tran, mock_ga_t *g_a, mock_ga_t *g_b);
void          Mock_Mask_sync(int first, int last);
void          Mock_Matmul_patch_2d(char transa, char transb, void* alpha, void *beta, mock_ga_t *g_a, int ailo, int aihi, int ajlo, int ajhi, mock_ga_t *g_b, int bilo, int bihi, int bjlo, int bjhi, mock_ga_t *g_c, int cilo, int cihi, int cjlo, int cjhi);
void          Mock_Matmul_patch(char transa, char transb, void* alpha, void *beta, mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[], mock_ga_t *g_c, int clo[], int chi[]) ;
void          Mock_Median(mock_ga_t *g_a, mock_ga_t *g_b, mock_ga_t *g_c, int g_m);
void          Mock_Median_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi, mock_ga_t *g_c, int *clo, int *chi, int g_m, int *mlo, int *mhi);
size_t        Mock_Memory_avail(void);
int           Mock_Memory_limited(void);
void          Mock_Merge_distr_patch(mock_ga_t *g_a, int alo[], int ahi[], mock_ga_t *g_b, int blo[], int bhi[]);
void          Mock_Merge_mirrored(mock_ga_t *g_a);
void          Mock_NbAcc(mock_ga_t *g_a,int lo[], int hi[],void* buf,int ld[],void* alpha, ga_nbhdl_t* nbhandle);
void          Mock_Nbget_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize,void *buf, int *ld, ga_nbhdl_t *nbhandle);
void          Mock_NbGet_ghost_dir(mock_ga_t *g_a, int mask[], ga_nbhdl_t* handle);
void          Mock_NbGet(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle);
void          Mock_Nblock(mock_ga_t *g_a, int *nblock);
void          Mock_Nbput_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld, ga_nbhdl_t *nbhandle);
void          Mock_NbPut(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[], ga_nbhdl_t* nbhandle);
int           Mock_NbTest(ga_nbhdl_t* nbhandle);
void          Mock_NbWait(ga_nbhdl_t* nbhandle);
int           Mock_Ndim(mock_ga_t *g_a);
int           Mock_Nnodes(void);
int           Mock_Nodeid(void);
void          Mock_Norm1(mock_ga_t *g_a, double *nm);
void          Mock_Norm_infinity(mock_ga_t *g_a, double *nm);
void          Mock_Periodic_acc(mock_ga_t *g_a, int lo[], int hi[],void* buf,int ld[],void* alpha);
void          Mock_Periodic_get(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[]); 
void          Mock_Periodic_put(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[]); 
int           Mock_Pgroup_absolute_id(int pgroup, int pid);
void          Mock_Pgroup_brdcst(int grp, void *buf, int lenbuf, int root);
void          Mock_Pgroup_cgop(int grp, SingleComplex x[], int n, char *op);
int           Mock_Pgroup_create(int *list, int count);
int           Mock_Pgroup_destroy(int grp);
void          Mock_Pgroup_dgop(int grp, double x[], int n, char *op);
void          Mock_Pgroup_fgop(int grp, float x[], int n, char *op);   
int           Mock_Pgroup_get_default(void);
int           Mock_Pgroup_get_mirror(void);
int           Mock_Pgroup_get_world(void);
void          Mock_Pgroup_igop(int grp, int x[], int n, char *op);
void          Mock_Pgroup_lgop(int grp, long x[], int n, char *op);
void          Mock_Pgroup_llgop(int grp, long long x[], int n, char *op);
int           Mock_Pgroup_nnodes(int grp_id);
int           Mock_Pgroup_nodeid(int grp_id);
void          Mock_Pgroup_set_default(int p_handle);
int           Mock_Pgroup_split(int grp_id, int num_group);
int           Mock_Pgroup_split_irreg(int grp_id, int color);
void          Mock_Pgroup_sync(int grp_id);
void          Mock_Pgroup_zgop(int grp, DoubleComplex x[], int n, char *op);
void          Mock_Print_distribution(mock_ga_t *g_a); 
void          Mock_Print_file(FILE *file, mock_ga_t *g_a);
void          Mock_Print(mock_ga_t *g_a);
void          Mock_Print_patch_2d(mock_ga_t *g_a,int ilo,int ihi,int jlo,int jhi,int pretty);
void          Mock_Print_patch(mock_ga_t *g_a, int lo[], int hi[], int pretty);
void          Mock_Print_stats(void);
void          Mock_Proc_topology(mock_ga_t *g_a, int proc, int coord[]);
void          Mock_Put_field(mock_ga_t *g_a, int *lo, int *hi, int foff, int fsize, void *buf, int *ld);
void          Mock_Put(mock_ga_t *g_a, int lo[], int hi[], void* buf, int ld[]); 
void          Mock_Randomize(mock_ga_t *g_a, void *value);
long          Mock_Read_inc(mock_ga_t *g_a, int subscript[], long inc);
void          Mock_Recip(mock_ga_t *g_a);
void          Mock_Recip_patch(mock_ga_t *g_a,int *lo, int *hi);
void          Mock_Register_stack_memory(void * (*ext_alloc)(size_t, int, char *), void (*ext_free)(void *));
int           Mock_Register_type(size_t bytes);
void          Mock_Release_block_grid(mock_ga_t *g_a, int index[]);
void          Mock_Release_block(mock_ga_t *g_a, int idx);
void          Mock_Release_block_segment(mock_ga_t *g_a, int idx);
void          Mock_Release_ghost_element(mock_ga_t *g_a, int index[]);
void          Mock_Release_ghosts(mock_ga_t *g_a);
void          Mock_Release(mock_ga_t *g_a, int lo[], int hi[]);
void          Mock_Release_update_block_grid(mock_ga_t *g_a, int index[]);
void          Mock_Release_update_block(mock_ga_t *g_a, int idx);
void          Mock_Release_update_block_segment(mock_ga_t *g_a, int idx);
void          Mock_Release_update_ghost_element(mock_ga_t *g_a, int index[]);
void          Mock_Release_update_ghosts(mock_ga_t *g_a);
void          Mock_Release_update(mock_ga_t *g_a, int lo[], int hi[]);
void          Mock_Scale_cols(mock_ga_t *g_a, int g_v);
void          Mock_Scale(mock_ga_t *g_a, void *value); 
void          Mock_Scale_patch(mock_ga_t *g_a, int lo[], int hi[], void *alpha);
void          Mock_Scale_rows(mock_ga_t *g_a, int g_v);
void          Mock_Scan_add(mock_ga_t *g_a, mock_ga_t *g_b, int g_sbit, int lo, int hi, int excl);
void          Mock_Scan_copy(mock_ga_t *g_a, mock_ga_t *g_b, int g_sbit, int lo, int hi);
void          Mock_Scatter_acc_flat(mock_ga_t *g_a, void *v, int subsArray[], int n, void *alpha);
void          Mock_Scatter_acc(mock_ga_t *g_a, void *v, int* subsArray[], int n, void *alpha);
void          Mock_Scatter_flat(mock_ga_t *g_a, void *v, int subsArray[], int n);
void          Mock_Scatter(mock_ga_t *g_a, void *v, int* subsArray[], int n);
void          Mock_Select_elem(mock_ga_t *g_a, char* op, void* val, int *index);
void          Mock_Set_array_name(mock_ga_t *g_a, char *name);
void          Mock_Set_block_cyclic(mock_ga_t *g_a, int dims[]);
void          Mock_Set_block_cyclic_proc_grid(mock_ga_t *g_a, int block[], int proc_grid[]);
void          Mock_Set_chunk(mock_ga_t *g_a, int chunk[]);
void          Mock_Set_data(mock_ga_t *g_a, int ndim, int dims[], int type);
void          Mock_Set_debug(int flag);
void          Mock_Set_diagonal(mock_ga_t *g_a, int g_v);
void          Mock_Set_ghost_corner_flag(mock_ga_t *g_a, int flag);
void          Mock_Set_ghosts(mock_ga_t *g_a, int width[]);
void          Mock_Set_irreg_distr(mock_ga_t *g_a, int map[], int block[]);
void          Mock_Set_irreg_flag(mock_ga_t *g_a, int flag);
void          Mock_Set_memory_limit(size_t limit);
void          Mock_Set_pgroup(mock_ga_t *g_a, int p_handle);
void          Mock_Set_restricted(mock_ga_t *g_a, int list[], int size);
void          Mock_Set_restricted_range(mock_ga_t *g_a, int lo_proc, int hi_proc);
void          Mock_Sgemm(char ta, char tb, int m, int n, int k, float alpha, mock_ga_t *g_a, mock_ga_t *g_b, float beta, mock_ga_t *g_c );
void          Mock_Shift_diagonal(mock_ga_t *g_a, void *c);
int           Mock_Solve(mock_ga_t *g_a, mock_ga_t *g_b);
int           Mock_Spd_invert(mock_ga_t *g_a);
void          Mock_Step_bound_info(int g_xx, int g_vv, int g_xxll, int g_xxuu, void *boundmin, void *wolfemin, void *boundmax);
void          Mock_Step_bound_info_patch(int g_xx, int *xxlo, int *xxhi, int g_vv, int *vvlo, int *vvhi, int g_xxll, int *xxlllo, int *xxllhi, int g_xxuu, int *xxuulo, int *xxuuhi, void *boundmin, void *wolfemin, void *boundmax);
void          Mock_Step_max(mock_ga_t *g_a, mock_ga_t *g_b, void *step);
void          Mock_Step_max_patch(mock_ga_t *g_a, int *alo, int *ahi, mock_ga_t *g_b, int *blo, int *bhi, void *step);
void          Mock_Strided_acc(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[], void *alpha); 
void          Mock_Strided_get(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[]); 
void          Mock_Strided_put(mock_ga_t *g_a, int lo[], int hi[], int skip[], void* buf, int ld[]); 
void          Mock_Summarize(int verbose);
void          Mock_Symmetrize(mock_ga_t *g_a);
void          Mock_Sync(void);
void          Mock_Terminate(void);
int           Mock_Total_blocks(mock_ga_t *g_a);   
void          Mock_Transpose(mock_ga_t *g_a, mock_ga_t *g_b);
void          Mock_Unlock(int mutex);
int           Mock_Update_ghost_dir(mock_ga_t *g_a, int dimension, int idir, int flag);
void          Mock_Update_ghosts(mock_ga_t *g_a);
int           Mock_Uses_fapi(void);
int           Mock_Uses_ma(void);
int           Mock_Uses_proc_grid(mock_ga_t *g_a);
int           Mock_Valid_handle(mock_ga_t *g_a);
int           Mock_Verify_handle(mock_ga_t *g_a);
double        Mock_Wtime(void);
DoubleComplex Mock_Zdot(mock_ga_t *g_a, mock_ga_t *g_b); 
DoubleComplex Mock_Zdot_patch(mock_ga_t *g_a, char t_a, int alo[], int ahi[], mock_ga_t *g_b, char t_b, int blo[], int bhi[]);
void          Mock_Zero_diagonal(mock_ga_t *g_a);
void          Mock_Zero(mock_ga_t *g_a); 
void          Mock_Zero_patch(mock_ga_t *g_a, int lo[], int hi[]);
void          Mock_Zgemm(char ta, char tb, int m, int n, int k, DoubleComplex alpha, mock_ga_t *g_a, mock_ga_t *g_b, DoubleComplex beta, mock_ga_t *g_c );
void          Mock_Zgop(DoubleComplex x[], int n, char *op);

#endif /* GLOBAL_TESTING_UNIT_TEST_MOCK_H */
