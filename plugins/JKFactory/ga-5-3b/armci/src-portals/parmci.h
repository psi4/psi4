#include "armci.h"

int PARMCI_AccV (int op, void *scale, armci_giov_t * darr, int len, int proc);

void PARMCI_Barrier ();

int PARMCI_AccS (int optype, void *scale, void *src_ptr, int *src_stride_arr,
		 void *dst_ptr, int *dst_stride_arr, int *count,
		 int stride_levels, int proc);

void PARMCI_Finalize ();

int PARMCI_NbPut (void *src, void *dst, int bytes, int proc,
		  armci_hdl_t * nb_handle);

int PARMCI_GetValueInt (void *src, int proc);

int PARMCI_Put_flag (void *src, void *dst, int bytes, int *f, int v,
		     int proc);

int PARMCI_NbGetS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		   int *dst_stride_arr, int *count, int stride_levels,
		   int proc, armci_hdl_t * nb_handle);

void *PARMCI_Malloc_local (armci_size_t bytes);

int PARMCI_Free_local (void *ptr);

int PARMCI_Get (void *src, void *dst, int bytes, int proc);

int PARMCI_Put (void *src, void *dst, int bytes, int proc);

int PARMCI_Destroy_mutexes ();

int PARMCI_GetS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		 int *dst_stride_arr, int *count, int stride_levels,
		 int proc);

int PARMCI_NbAccV (int op, void *scale, armci_giov_t * darr, int len,
		   int proc, armci_hdl_t * nb_handle);

float PARMCI_GetValueFloat (void *src, int proc);

int PARMCI_Malloc (void **ptr_arr, armci_size_t bytes);

int PARMCI_NbAccS (int optype, void *scale, void *src_ptr,
		   int *src_stride_arr, void *dst_ptr, int *dst_stride_arr,
		   int *count, int stride_levels, int proc,
		   armci_hdl_t * nb_handle);

int PARMCI_PutS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		 int *dst_stride_arr, int *count, int stride_levels,
		 int proc);

void *PARMCI_Memat (armci_meminfo_t * meminfo, long offset);

int PARMCI_PutV (armci_giov_t * darr, int len, int proc);

int PARMCI_Free (void *ptr);

int PARMCI_Init_args (int *argc, char ***argv);

int PARMCI_PutValueInt (int src, void *dst, int proc);

void PARMCI_Memget (size_t bytes, armci_meminfo_t * meminfo, int memflg);

void PARMCI_AllFence ();

int PARMCI_NbPutV (armci_giov_t * darr, int len, int proc,
		   armci_hdl_t * nb_handle);

int PARMCI_PutValueDouble (double src, void *dst, int proc);

int PARMCI_GetV (armci_giov_t * darr, int len, int proc);

int PARMCI_Test (armci_hdl_t * nb_handle);

void PARMCI_Unlock (int mutex, int proc);

void PARMCI_Fence (int proc);

int PARMCI_Create_mutexes (int num);

int PARMCI_PutS_flag (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		      int *dst_stride_arr, int *count, int stride_levels,
		      int *flag, int val, int proc);

int PARMCI_WaitProc (int proc);

void PARMCI_Lock (int mutex, int proc);

double PARMCI_GetValueDouble (void *src, int proc);

int PARMCI_NbGetV (armci_giov_t * darr, int len, int proc,
		   armci_hdl_t * nb_handle);

int PARMCI_Rmw (int op, int *ploc, int *prem, int extra, int proc);

int PARMCI_Init ();

int PARMCI_WaitAll ();

int PARMCI_NbGet (void *src, void *dst, int bytes, int proc,
		  armci_hdl_t * nb_handle);

int PARMCI_PutValueFloat (float src, void *dst, int proc);

int PARMCI_NbPutS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		   int *dst_stride_arr, int *count, int stride_levels,
		   int proc, armci_hdl_t * nb_handle);

int PARMCI_PutS_flag_dir (void *src_ptr, int *src_stride_arr, void *dst_ptr,
			  int *dst_stride_arr, int *count, int stride_levels,
			  int *flag, int val, int proc);

int PARMCI_PutValueLong (long src, void *dst, int proc);

int PARMCI_Wait (armci_hdl_t * nb_handle);

long PARMCI_GetValueLong (void *src, int proc);
