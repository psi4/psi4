#if HAVE_CONFIG_H
#   include "config.h"
#endif


#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#include "armci.h"
#include "parmci.h"

int
ARMCI_AccV (int op, void *scale, armci_giov_t * darr, int len, int proc)
{
  int rval;
  rval = PARMCI_AccV (op, scale, darr, len, proc);
  return rval;
}

void
ARMCI_Barrier ()
{
  PARMCI_Barrier ();
}

int
ARMCI_AccS (int optype, void *scale, void *src_ptr, int *src_stride_arr,
	    void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels,
	    int proc)
{
  int rval;
  rval =
    PARMCI_AccS (optype, scale, src_ptr, src_stride_arr, dst_ptr,
		 dst_stride_arr, count, stride_levels, proc);
  return rval;
}

void
ARMCI_Finalize ()
{
  PARMCI_Finalize ();
}

int
ARMCI_NbPut (void *src, void *dst, int bytes, int proc,
	     armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_NbPut (src, dst, bytes, proc, nb_handle);
  return rval;
}

int
ARMCI_GetValueInt (void *src, int proc)
{
  int rval;
  rval = PARMCI_GetValueInt (src, proc);
  return rval;
}

int
ARMCI_Put_flag (void *src, void *dst, int bytes, int *f, int v, int proc)
{
  int rval;
  rval = PARMCI_Put_flag (src, dst, bytes, f, v, proc);
  return rval;
}

int
ARMCI_NbGetS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
	      int *dst_stride_arr, int *count, int stride_levels, int proc,
	      armci_hdl_t * nb_handle)
{
  int rval;
  rval =
    PARMCI_NbGetS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		   stride_levels, proc, nb_handle);
  return rval;
}

void *
ARMCI_Malloc_local (armci_size_t bytes)
{
  void *rval;
  rval = PARMCI_Malloc_local (bytes);
  return rval;
}

int
ARMCI_Free_local (void *ptr)
{
  int rval;
  rval = PARMCI_Free_local (ptr);
  return rval;
}

int
ARMCI_Get (void *src, void *dst, int bytes, int proc)
{
  int rval;
  rval = PARMCI_Get (src, dst, bytes, proc);
  return rval;
}

int
ARMCI_Put (void *src, void *dst, int bytes, int proc)
{
  int rval;
  rval = PARMCI_Put (src, dst, bytes, proc);
  return rval;
}

int
ARMCI_Destroy_mutexes ()
{
  int rval;
  rval = PARMCI_Destroy_mutexes ();
  return rval;
}

int
ARMCI_GetS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
	    int *dst_stride_arr, int *count, int stride_levels, int proc)
{
  int rval;
  rval =
    PARMCI_GetS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		 stride_levels, proc);
  return rval;
}

int
ARMCI_NbAccV (int op, void *scale, armci_giov_t * darr, int len, int proc,
	      armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_NbAccV (op, scale, darr, len, proc, nb_handle);
  return rval;
}

float
ARMCI_GetValueFloat (void *src, int proc)
{
  float rval;
  rval = PARMCI_GetValueFloat (src, proc);
  return rval;
}

int
ARMCI_Malloc (void **ptr_arr, armci_size_t bytes)
{
  int rval;
  rval = PARMCI_Malloc (ptr_arr, bytes);
  return rval;
}

int
ARMCI_NbAccS (int optype, void *scale, void *src_ptr, int *src_stride_arr,
	      void *dst_ptr, int *dst_stride_arr, int *count,
	      int stride_levels, int proc, armci_hdl_t * nb_handle)
{
  int rval;
  rval =
    PARMCI_NbAccS (optype, scale, src_ptr, src_stride_arr, dst_ptr,
		   dst_stride_arr, count, stride_levels, proc, nb_handle);
  return rval;
}

int
ARMCI_PutS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
	    int *dst_stride_arr, int *count, int stride_levels, int proc)
{
  int rval;
  rval =
    PARMCI_PutS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		 stride_levels, proc);
  return rval;
}

void *
ARMCI_Memat (armci_meminfo_t * meminfo, long offset)
{
  void *rval;
  rval = PARMCI_Memat (meminfo, offset);
  return rval;
}

int
ARMCI_PutV (armci_giov_t * darr, int len, int proc)
{
  int rval;
  rval = PARMCI_PutV (darr, len, proc);
  return rval;
}

int
ARMCI_Free (void *ptr)
{
  int rval;
  rval = PARMCI_Free (ptr);
  return rval;
}

int
ARMCI_Init_args (int *argc, char ***argv)
{
  int rval;
  rval = PARMCI_Init_args (argc, argv);
  return rval;
}

int
ARMCI_PutValueInt (int src, void *dst, int proc)
{
  int rval;
  rval = PARMCI_PutValueInt (src, dst, proc);
  return rval;
}

void
ARMCI_Memget (size_t bytes, armci_meminfo_t * meminfo, int memflg)
{
  PARMCI_Memget (bytes, meminfo, memflg);
}

void
ARMCI_AllFence ()
{
  PARMCI_AllFence ();
}

int
ARMCI_NbPutV (armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_NbPutV (darr, len, proc, nb_handle);
  return rval;
}

int
ARMCI_PutValueDouble (double src, void *dst, int proc)
{
  int rval;
  rval = PARMCI_PutValueDouble (src, dst, proc);
  return rval;
}

int
ARMCI_GetV (armci_giov_t * darr, int len, int proc)
{
  int rval;
  rval = PARMCI_GetV (darr, len, proc);
  return rval;
}

int
ARMCI_Test (armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_Test (nb_handle);
  return rval;
}

void
ARMCI_Unlock (int mutex, int proc)
{
  PARMCI_Unlock (mutex, proc);
}

void
ARMCI_Fence (int proc)
{
  PARMCI_Fence (proc);
}

int
ARMCI_Create_mutexes (int num)
{
  int rval;
  rval = PARMCI_Create_mutexes (num);
  return rval;
}

int
ARMCI_PutS_flag (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		 int *dst_stride_arr, int *count, int stride_levels,
		 int *flag, int val, int proc)
{
  int rval;
  rval =
    PARMCI_PutS_flag (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		      stride_levels, flag, val, proc);
  return rval;
}

int
ARMCI_WaitProc (int proc)
{
  int rval;
  rval = PARMCI_WaitProc (proc);
  return rval;
}

void
ARMCI_Lock (int mutex, int proc)
{
  PARMCI_Lock (mutex, proc);
}

double
ARMCI_GetValueDouble (void *src, int proc)
{
  double rval;
  rval = PARMCI_GetValueDouble (src, proc);
  return rval;
}

int
ARMCI_NbGetV (armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_NbGetV (darr, len, proc, nb_handle);
  return rval;
}

int
ARMCI_Rmw (int op, int *ploc, int *prem, int extra, int proc)
{
  int rval;
  rval = PARMCI_Rmw (op, ploc, prem, extra, proc);
  return rval;
}

int
ARMCI_Init ()
{
  int rval;
  rval = PARMCI_Init ();
  return rval;
}

int
ARMCI_WaitAll ()
{
  int rval;
  rval = PARMCI_WaitAll ();
  return rval;
}

int
ARMCI_NbGet (void *src, void *dst, int bytes, int proc,
	     armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_NbGet (src, dst, bytes, proc, nb_handle);
  return rval;
}

int
ARMCI_PutValueFloat (float src, void *dst, int proc)
{
  int rval;
  rval = PARMCI_PutValueFloat (src, dst, proc);
  return rval;
}

int
ARMCI_NbPutS (void *src_ptr, int *src_stride_arr, void *dst_ptr,
	      int *dst_stride_arr, int *count, int stride_levels, int proc,
	      armci_hdl_t * nb_handle)
{
  int rval;
  rval =
    PARMCI_NbPutS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		   stride_levels, proc, nb_handle);
  return rval;
}

int
ARMCI_PutS_flag_dir (void *src_ptr, int *src_stride_arr, void *dst_ptr,
		     int *dst_stride_arr, int *count, int stride_levels,
		     int *flag, int val, int proc)
{
  int rval;
  rval =
    PARMCI_PutS_flag_dir (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
			  count, stride_levels, flag, val, proc);
  return rval;
}

int
ARMCI_PutValueLong (long src, void *dst, int proc)
{
  int rval;
  rval = PARMCI_PutValueLong (src, dst, proc);
  return rval;
}

int
ARMCI_Wait (armci_hdl_t * nb_handle)
{
  int rval;
  rval = PARMCI_Wait (nb_handle);
  return rval;
}

long
ARMCI_GetValueLong (void *src, int proc)
{
  long rval;
  rval = PARMCI_GetValueLong (src, proc);
  return rval;
}
