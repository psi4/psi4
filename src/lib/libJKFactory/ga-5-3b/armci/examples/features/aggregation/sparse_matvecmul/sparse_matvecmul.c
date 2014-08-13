#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id$ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include "armci.h"
#include "message.h"

#define DIM1 5
#define DIM2 3
#ifdef __sun
/* Solaris has shared memory shortages in the default system configuration */
# define DIM3 6
# define DIM4 5
# define DIM5 4
#elif defined(__alpha__)
# define DIM3 8
# define DIM4 5
# define DIM5 6
#else
# define DIM3 8
# define DIM4 9
# define DIM5 7
#endif
#define DIM6 3
#define DIM7 2


#define OFF 1
#define EDIM1 (DIM1+OFF)
#define EDIM2 (DIM2+OFF)
#define EDIM3 (DIM3+OFF)
#define EDIM4 (DIM4+OFF)
#define EDIM5 (DIM5+OFF)
#define EDIM6 (DIM6+OFF)
#define EDIM7 (DIM7+OFF)

#define DIMS 4
#define MAXDIMS 7
#define MAX_DIM_VAL 50 
#define LOOP 200

#define BASE 100.
#define MAXPROC 128
#define TIMES 100

#ifdef CRAY
# define ELEMS 800
#else
# define ELEMS 200
#endif

/***************************** macros ************************/
#define COPY(src, dst, bytes) memcpy((dst),(src),(bytes))
#define ARMCI_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define ARMCI_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ARMCI_ABS(a) (((a) <0) ? -(a) : (a))

#define ROW      65536
#define COL      ROW   /* square matrices only for the time being */

/***************************** global data *******************/
int me, nproc;
short int fortran_indexing=0;
static int proc_row_list[MAXPROC];/*no of rows owned by each process - accumulated*/
static int proc_nz_list[MAXPROC]; /*no of non-zeros owned by each process */

#ifdef MSG_COMMS_PVM
void pvm_init(int argc, char *argv[])
{
    int mytid, mygid, ctid[MAXPROC];
    int np, i;

    mytid = pvm_mytid();
    if((argc != 2) && (argc != 1)) goto usage;
    if(argc == 1) np = 1;
    if(argc == 2)
        if((np = atoi(argv[1])) < 1) goto usage;
    if(np > MAXPROC) goto usage;

    mygid = pvm_joingroup(MPGROUP);

    if(np > 1)
        if (mygid == 0) 
            i = pvm_spawn(argv[0], argv+1, 0, "", np-1, ctid);

    while(pvm_gsize(MPGROUP) < np) sleep(1);

    /* sync */
    pvm_barrier(MPGROUP, np);
    
    printf("PVM initialization done!\n");
    
    return;

usage:
    fprintf(stderr, "usage: %s <nproc>\n", argv[0]);
    pvm_exit();
    exit(-1);
}
#endif
          
void create_array(void *a[], int elem_size, int ndim, int dims[])
{
     int bytes=elem_size, i, rc;

     assert(ndim<=MAXDIMS);
     for(i=0;i<ndim;i++)bytes*=dims[i];

     rc = ARMCI_Malloc(a, bytes);
     assert(rc==0);
     
     assert(a[me]);
     
}

void destroy_array(void *ptr[])
{
    armci_msg_barrier();

    assert(!ARMCI_Free(ptr[me]));
}

static void verify_list(int *proc_row_list) {
  int i;
  printf("\nVERIFY: %d: No of rows = %d\n\n", 0, proc_row_list[0]);
  for(i=1; i<nproc; i++)
    printf("\nVERIFY: %d: No of rows = %d\n\n", i, proc_row_list[i]-proc_row_list[i-1]);    
  fflush(stdout);
}

static void load_balance(int n, int non_zero, int *row_ind_tmp) {

  int proc_id, i, local_nz, local_nz_acc, A, B;

  local_nz = local_nz_acc = non_zero/nproc;

  /* number of rows owned by each process is stored in proc_row_list. This 
     is supposed to be well load balanced, so that each process has almost 
     same number of non-zero elements */
  proc_id = 0;
  if(me==0) printf("local_nz = %d\n", local_nz);
  for(i=0; i<n; i++) { /* as # of entries in row_ind_tmp = n+1 */
    if(row_ind_tmp[i] < local_nz_acc && row_ind_tmp[i+1] >= local_nz_acc) {
      proc_row_list[proc_id++] = i+1;
      local_nz_acc = local_nz*(proc_id+1);
      if(proc_id == nproc-1) local_nz_acc = non_zero;
      if(me==0 && proc_id<nproc) printf("local_nz = %d\n", local_nz_acc);
    }
  }

  proc_row_list[nproc-1] = n;

  for(i=0; i<nproc; i++) {
    A = (i==0) ? 0: proc_row_list[i-1];/* # of entries in row_ind_tmp is n+1*/ 
    B = proc_row_list[i];
    proc_nz_list[i] = row_ind_tmp[B]-row_ind_tmp[A];
  }
  
  if(proc_id != nproc) 
    ARMCI_Error("Error while preparing Process Row list", proc_id-1);

#if 1
  if(me==0) verify_list(proc_row_list);
#endif

}

static int sparse_initialize(int *n, int *non_zero, int **row_ind, 
                 int **col_ind, double **values, double **vec,
                 double **svec) {
  
  int i, j, rc, max, *row_ind_tmp=NULL, *tmp_indices=NULL;
  double *tmp_values=NULL;
  unsigned long len;
  FILE *fp=NULL;

  /* Broadcast order of matrix */
  if(me==0) {
    if((fp=fopen("Sparse-MPI/av41092.rua.data", "r")) == NULL)
      ARMCI_Error("Error: Input file not found", me);
    fortran_indexing = 1; /* This is 1 for Harwell-Boeing format matrices */
    fscanf(fp, "%d", n);
    if(*n%nproc) 
      ARMCI_Error("# of rows is not divisible by # of processors", nproc);
    if(*n > ROW) 
      ARMCI_Error("order is greater than defined variable ROW", ROW);
  }
  len = sizeof(int);
  armci_msg_brdcst(n, len, 0);  

  /* Broad cast number of non_zeros */
  if(me==0) fscanf(fp, "%d", non_zero);
  armci_msg_brdcst(non_zero, len, 0); 

  /* Broadcast row indices */
  len = (*n+1)*sizeof(int);
  row_ind_tmp = (int *)malloc(len);
  if(me==0)for(i=0; i<*n+1; i++) {
    fscanf(fp, "%d", &row_ind_tmp[i]);
    if(fortran_indexing) --row_ind_tmp[i];
  }
  armci_msg_brdcst(row_ind_tmp, len, 0);  
  
  load_balance(*n, *non_zero, row_ind_tmp);
  
  /* find how much temporary storage is needed at the maximum */
  if(me==0) {
    for(max=-1,j=0;j<nproc;j++) if(max<proc_nz_list[j]) max=proc_nz_list[j];
    if(max<0) ARMCI_Error(" max cannot be negative", max);
  }
  
  /* Broadcast the maximum number of elements */
  len = sizeof(int);
  armci_msg_brdcst(&max, len, 0); 

  /* create the Sparse MAtrix Array */
  if(me==0) printf("  Creating ValueArray (CompressedSparseMatrix) ...\n\n");
  create_array((void**)col_ind, sizeof(int), 1, &max);
   
  /* create the column subscript array */
  if(me==0) printf("  Creating Column Subscript Array ... \n\n");
  create_array((void**)values, sizeof(double), 1, &max);

  /* create the x-vector and the solution vector */
  if(me==0) printf("  Creating Vectors ... \n\n");
  create_array((void**)vec,  sizeof(double),1, &max);
  create_array((void**)svec, sizeof(double),1, &max);
  armci_msg_barrier();

  
  /* Process 0 distributes the column indices and non_zero values to 
     respective processors*/
  if(me == 0) {
    tmp_indices = (int *)malloc(max*sizeof(int));
    tmp_values  = (double *)malloc(max*sizeof(double));
    
    for(j=0; j<nproc; j++) {
      for(i=0; i<proc_nz_list[j]; i++) {
    fscanf(fp, "%d", &tmp_indices[i]); 
    if(fortran_indexing) --tmp_indices[i];
      }
      /* rc = fread(tmp_indices, sizeof(int), proc_nz_list[j], fp); */
      if((rc=ARMCI_Put(tmp_indices, col_ind[j], proc_nz_list[j]*sizeof(int), j)))
    ARMCI_Error("armci_nbput failed\n",rc);
    }
    for(j=0; j<nproc; j++) {
      for(i=0; i<proc_nz_list[j]; i++) fscanf(fp, "%lf", &tmp_values[i]);
      if((rc=ARMCI_Put(tmp_values, values[j], proc_nz_list[j]*sizeof(double), j)))
    ARMCI_Error("armci_nbput failed\n",rc);
    }
  }
  ARMCI_AllFence(); armci_msg_barrier();ARMCI_AllFence();

  /* initializing x-vector */
  if(me==0) for(i=0;i<proc_nz_list[me]; i++) vec[me][i] = (i+1);
  else for(i=0;i<proc_nz_list[me];i++) vec[me][i]=me*proc_nz_list[me-1]+(i+1);

#if 0
  if(me==0) {
    printf("max = %d\n", max);
    for(i=0; i<max; i++)  printf("%.1f ", values[me][i]);
    printf("\n");
  }
#endif

  *row_ind = row_ind_tmp;
  if(me==0) {
    free(tmp_indices);
    free(tmp_values);
    fclose(fp);
  }
  return 0;
}

static int compare(const void *p1, const void *p2) {
  int i = *((int *)p1);
  int j = *((int *)p2);
  
  if (i > j) return (1);
  if (i < j) return (-1);
  return (0);
}

static int count = -1;
static armci_hdl_t gHandle[MAXPROC];
static int prev_proc = -1;

static void get_data(int n, int start, int end, double *vec_local, 
            double **vec) {
  int i, j, rc, bytes, offset;
  int proc_start, proc_end, idx_start, idx_end;

  proc_start = proc_end = -1;
  for(i=0; i<nproc; i++) {
    if(proc_start<0 && proc_row_list[i]>start) proc_start = i;
    if(proc_end<0 && proc_row_list[i]>end) proc_end = i;    
  }
  if(proc_start<0 || proc_end<0) ARMCI_Error("Invalid Process Ids", -1);

  for(i=proc_start; i<=proc_end; i++) {
    if(i==proc_start) idx_start = start;
    else { if(i==0) idx_start=0; else idx_start = proc_row_list[i-1];}
    if(i==proc_end) idx_end = end;
    else idx_end = proc_row_list[i]-1;
   
    if(i!=prev_proc) {
      ++count;   prev_proc = i;  
      ARMCI_INIT_HANDLE(&gHandle[count]);
      ARMCI_SET_AGGREGATE_HANDLE(&gHandle[count]);
    }
    
    if(i==0) offset=0; else offset = proc_row_list[i-1];
    if(i==me) { /* local */
      for(j=idx_start; j<=idx_end; j++) vec_local[j] = vec[me][j-offset];
    }
    else {     /* remote */
      bytes = (idx_end-idx_start+1)*sizeof(double);
      vec_local[idx_start] = -1;
#if 0
      if((rc=ARMCI_Get(&vec[i][idx_start-offset], &vec_local[idx_start],
               bytes, i)))
#else
      if((rc=ARMCI_NbGet(&vec[i][idx_start-offset], &vec_local[idx_start],
               bytes, i, &gHandle[count])))
#endif
    ARMCI_Error("armci_nbget failed\n",rc);
    }
  }
}

static void sparse_multiply(int n, int non_zero, int *row_ind, int **col_ind, 
            double **values, double **vec, double **svec) {
  
  int i, j, k, num_elements, offset, *tmp_indices;
  double start_time, comm_time, comp_time, v, vec_local[COL];
  int start, end, prev, nrows, idx;

#if 0
  /* ---- Sequential Case ----  */
   for(i=0; i<n; i++) {
     svec[me][i] = 0;
     for(k=row_ind[i]; k<row_ind[i+1]; k++) {
       j = col_ind[me][k];
       v = values[me][k];
       svec[me][i] += v*vec[me][j];
       printf("%.1f %.1f\n", v, vec[me][j]);
     }
   }
   for(i=0; i<n; i++) printf("%.1f ", svec[me][i]);
   printf("\n");
#else

  num_elements = proc_nz_list[me];
  printf("num_elements = %d\n", num_elements);
  tmp_indices = (int *)malloc(num_elements*sizeof(int));
  for(i=0; i<num_elements; i++) tmp_indices[i] = col_ind[me][i];
  qsort(tmp_indices, num_elements, sizeof(int), compare);

  start_time = armci_timer();

  /* get the required portion of vector you need to local array */
  start = prev = tmp_indices[0];
  for(i=1; i<num_elements; i++) {
    if(tmp_indices[i]>prev+1) {
      end = prev;
      get_data(n, start, end, vec_local, vec);
      start = prev = tmp_indices[i];
    }
    else prev = tmp_indices[i];
  }
  get_data(n, start, prev, vec_local, vec);

#if 1
  if(count>=0) for(i=0; i<=count; i++) ARMCI_Wait(&gHandle[i]);
#endif

  comm_time = armci_timer() - start_time;
  start_time = armci_timer();   

  /* Perform Matrix-Vector multiply and store the result in
     solution vector - "svec[]" */

  if(me==0) { nrows = proc_row_list[me]; offset = row_ind[0]; }
  else { 
    nrows = proc_row_list[me]-proc_row_list[me-1]; 
    offset = row_ind[proc_row_list[me-1]]; 
  }
  /* printf("%d: My total Work = %d\n", me, nrows); */

  for(i=0; i<nrows; i++) { /* loop over rows owned by me */
    svec[me][i] = 0;
    if(me==0) idx = i; else idx = proc_row_list[me-1] + i;
    for(k=row_ind[idx]; k<row_ind[idx+1]; k++) {
      j = col_ind[me][k-offset];
      v = values[me][k-offset];
      svec[me][i] += v*vec_local[j];
    }
  }
  comp_time = armci_timer()-start_time;
  printf("%d: %f + %f = %f  (count = %d)\n", me, comm_time, comp_time, 
     comm_time+comp_time, count+1);
#endif
}

static void gather_solution_vector(double **svec) {
#if 0
  double y[COL];
  if((rc=ARMCI_Get(&vec[i][idx_start-offset], &vec_local[idx_start],
           bytes, i)))
    ARMCI_Error("armci_nbget failed\n",rc);
#endif
}

static void test_sparse() {
  
    int *col_ind[MAXPROC];
    double *values[MAXPROC], *vec[MAXPROC], *svec[MAXPROC]/*, start_time*/;
    int n, non_zero, *row_ind;

    sparse_initialize(&n, &non_zero, &row_ind, col_ind, values, vec, svec);
    armci_msg_barrier();

    /*start_time = armci_timer();*/
    sparse_multiply(n, non_zero, row_ind, col_ind, values, vec, svec);
    /* printf("%d: Timetaken = %f\n", me, armci_timer()-start_time); */
    armci_msg_barrier();
    
    if(me==0) gather_solution_vector(svec);
    
    if(me==0){printf("O.K.\n"); fflush(stdout);}
    destroy_array((void **)vec);
}


int main(int argc, char* argv[])
{

    armci_msg_init(&argc, &argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();

/*    printf("nproc = %d, me = %d\n", nproc, me);*/
    
    if(nproc>MAXPROC && me==0)
       ARMCI_Error("Test works for up to %d processors\n",MAXPROC);

    if(me==0){
       printf("ARMCI test program (%d processes)\n",nproc); 
       fflush(stdout);
       sleep(1);
    }
    
    ARMCI_Init();

    if(me==0){
      printf("\n  Performing Sparse Matrix-Vector Multiplication ...\n\n");
      fflush(stdout);
    }
    test_sparse();
    
    ARMCI_AllFence();
    armci_msg_barrier();
    if(me==0){printf("\nSuccess!!\n"); fflush(stdout);}
    sleep(2);
    
    armci_msg_barrier();
    ARMCI_Finalize();
    armci_msg_finalize();
    return(0);
}
