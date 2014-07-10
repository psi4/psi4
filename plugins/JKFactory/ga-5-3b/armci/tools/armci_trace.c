#include <stdio.h>

#include <mpi.h>

#include "parmci.h"

#define TIME MPI_Wtime

static double ARMCI_AccV_t;
static double ARMCI_Barrier_t;
static double ARMCI_AccS_t;
static double ARMCI_Finalize_t;
static double ARMCI_NbPut_t;
static double ARMCI_GetValueInt_t;
static double ARMCI_Put_flag_t;
static double ARMCI_WaitAll_t;
static double ARMCI_Malloc_local_t;
static double ARMCI_Free_local_t;
static double ARMCI_Get_t;
static double ARMCI_PutValueFloat_t;
static double ARMCI_NbAccV_t;
static double ARMCI_GetValueFloat_t;
static double ARMCI_Malloc_t;
static double ARMCI_NbAccS_t;
static double ARMCI_PutS_t;
static double ARMCI_PutV_t;
static double ARMCI_Destroy_mutexes_t;
static double ARMCI_Free_t;
static double ARMCI_Init_args_t;
static double ARMCI_PutValueInt_t;
static double ARMCI_Memget_t;
static double ARMCI_AllFence_t;
static double ARMCI_NbPutV_t;
static double ARMCI_PutValueDouble_t;
static double ARMCI_GetV_t;
static double ARMCI_Test_t;
static double ARMCI_GetS_t;
static double ARMCI_Unlock_t;
static double ARMCI_Fence_t;
static double ARMCI_Create_mutexes_t;
static double ARMCI_PutS_flag_t;
static double ARMCI_WaitProc_t;
static double ARMCI_Lock_t;
static double ARMCI_GetValueDouble_t;
static double ARMCI_NbGetV_t;
static double ARMCI_Rmw_t;
static double ARMCI_Init_t;
static double ARMCI_NbGetS_t;
static double ARMCI_NbGet_t;
static double ARMCI_Put_t;
static double ARMCI_NbPutS_t;
static double ARMCI_PutS_flag_dir_t;
static double ARMCI_Wait_t;
static double ARMCI_GetValueLong_t;
static double ARMCI_PutValueLong_t;
static int ARMCI_AccV_c;
static int ARMCI_Barrier_c;
static int ARMCI_AccS_c;
static int ARMCI_Finalize_c;
static int ARMCI_NbPut_c;
static int ARMCI_GetValueInt_c;
static int ARMCI_Put_flag_c;
static int ARMCI_WaitAll_c;
static int ARMCI_Malloc_local_c;
static int ARMCI_Free_local_c;
static int ARMCI_Get_c;
static int ARMCI_PutValueFloat_c;
static int ARMCI_NbAccV_c;
static int ARMCI_GetValueFloat_c;
static int ARMCI_Malloc_c;
static int ARMCI_NbAccS_c;
static int ARMCI_PutS_c;
static int ARMCI_PutV_c;
static int ARMCI_Destroy_mutexes_c;
static int ARMCI_Free_c;
static int ARMCI_Init_args_c;
static int ARMCI_PutValueInt_c;
static int ARMCI_Memget_c;
static int ARMCI_AllFence_c;
static int ARMCI_NbPutV_c;
static int ARMCI_PutValueDouble_c;
static int ARMCI_GetV_c;
static int ARMCI_Test_c;
static int ARMCI_GetS_c;
static int ARMCI_Unlock_c;
static int ARMCI_Fence_c;
static int ARMCI_Create_mutexes_c;
static int ARMCI_PutS_flag_c;
static int ARMCI_WaitProc_c;
static int ARMCI_Lock_c;
static int ARMCI_GetValueDouble_c;
static int ARMCI_NbGetV_c;
static int ARMCI_Rmw_c;
static int ARMCI_Init_c;
static int ARMCI_NbGetS_c;
static int ARMCI_NbGet_c;
static int ARMCI_Put_c;
static int ARMCI_NbPutS_c;
static int ARMCI_PutS_flag_dir_c;
static int ARMCI_Wait_c;
static int ARMCI_GetValueLong_c;
static int ARMCI_PutValueLong_c;

static double t;
static int c;
static int me;

void ARMCI_Finalize()
{

    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    MPI_Reduce(&ARMCI_AccV_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_AccV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_AccV,%d,%lf\n", ARMCI_AccV_c, ARMCI_AccV_t);
    }

    MPI_Reduce(&ARMCI_Barrier_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Barrier_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Barrier,%d,%lf\n", ARMCI_Barrier_c, ARMCI_Barrier_t);
    }

    MPI_Reduce(&ARMCI_AccS_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_AccS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_AccS,%d,%lf\n", ARMCI_AccS_c, ARMCI_AccS_t);
    }

    MPI_Reduce(&ARMCI_NbPut_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbPut_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbPut,%d,%lf\n", ARMCI_NbPut_c, ARMCI_NbPut_t);
    }

    MPI_Reduce(&ARMCI_GetValueInt_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetValueInt_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetValueInt,%d,%lf\n", ARMCI_GetValueInt_c,
	       ARMCI_GetValueInt_t);
    }

    MPI_Reduce(&ARMCI_Put_flag_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Put_flag_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Put_flag,%d,%lf\n", ARMCI_Put_flag_c,
	       ARMCI_Put_flag_t);
    }

    MPI_Reduce(&ARMCI_NbGetS_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbGetS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbGetS,%d,%lf\n", ARMCI_NbGetS_c, ARMCI_NbGetS_t);
    }

    MPI_Reduce(&ARMCI_Malloc_local_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Malloc_local_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Malloc_local,%d,%lf\n", ARMCI_Malloc_local_c,
	       ARMCI_Malloc_local_t);
    }

    MPI_Reduce(&ARMCI_Free_local_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Free_local_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Free_local,%d,%lf\n", ARMCI_Free_local_c,
	       ARMCI_Free_local_t);
    }

    MPI_Reduce(&ARMCI_Get_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Get_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Get,%d,%lf\n", ARMCI_Get_c, ARMCI_Get_t);
    }

    MPI_Reduce(&ARMCI_Put_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Put_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Put,%d,%lf\n", ARMCI_Put_c, ARMCI_Put_t);
    }

    MPI_Reduce(&ARMCI_Destroy_mutexes_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Destroy_mutexes_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Destroy_mutexes,%d,%lf\n", ARMCI_Destroy_mutexes_c,
	       ARMCI_Destroy_mutexes_t);
    }

    MPI_Reduce(&ARMCI_GetS_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetS,%d,%lf\n", ARMCI_GetS_c, ARMCI_GetS_t);
    }

    MPI_Reduce(&ARMCI_NbAccV_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbAccV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbAccV,%d,%lf\n", ARMCI_NbAccV_c, ARMCI_NbAccV_t);
    }

    MPI_Reduce(&ARMCI_GetValueFloat_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetValueFloat_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetValueFloat,%d,%lf\n", ARMCI_GetValueFloat_c,
	       ARMCI_GetValueFloat_t);
    }

    MPI_Reduce(&ARMCI_Malloc_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Malloc_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Malloc,%d,%lf\n", ARMCI_Malloc_c, ARMCI_Malloc_t);
    }

    MPI_Reduce(&ARMCI_NbAccS_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbAccS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbAccS,%d,%lf\n", ARMCI_NbAccS_c, ARMCI_NbAccS_t);
    }

    MPI_Reduce(&ARMCI_PutS_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutS,%d,%lf\n", ARMCI_PutS_c, ARMCI_PutS_t);
    }

    MPI_Reduce(&ARMCI_PutV_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutV,%d,%lf\n", ARMCI_PutV_c, ARMCI_PutV_t);
    }

    MPI_Reduce(&ARMCI_Free_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Free_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Free,%d,%lf\n", ARMCI_Free_c, ARMCI_Free_t);
    }

    MPI_Reduce(&ARMCI_Init_args_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Init_args_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Init_args,%d,%lf\n", ARMCI_Init_args_c,
	       ARMCI_Init_args_t);
    }

    MPI_Reduce(&ARMCI_PutValueInt_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutValueInt_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutValueInt,%d,%lf\n", ARMCI_PutValueInt_c,
	       ARMCI_PutValueInt_t);
    }

    MPI_Reduce(&ARMCI_Memget_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Memget_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Memget,%d,%lf\n", ARMCI_Memget_c, ARMCI_Memget_t);
    }

    MPI_Reduce(&ARMCI_AllFence_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_AllFence_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_AllFence,%d,%lf\n", ARMCI_AllFence_c,
	       ARMCI_AllFence_t);
    }

    MPI_Reduce(&ARMCI_NbPutV_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbPutV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbPutV,%d,%lf\n", ARMCI_NbPutV_c, ARMCI_NbPutV_t);
    }

    MPI_Reduce(&ARMCI_PutValueDouble_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutValueDouble_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutValueDouble,%d,%lf\n", ARMCI_PutValueDouble_c,
	       ARMCI_PutValueDouble_t);
    }

    MPI_Reduce(&ARMCI_GetV_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetV,%d,%lf\n", ARMCI_GetV_c, ARMCI_GetV_t);
    }

    MPI_Reduce(&ARMCI_Test_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Test_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Test,%d,%lf\n", ARMCI_Test_c, ARMCI_Test_t);
    }

    MPI_Reduce(&ARMCI_Unlock_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Unlock_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Unlock,%d,%lf\n", ARMCI_Unlock_c, ARMCI_Unlock_t);
    }

    MPI_Reduce(&ARMCI_Fence_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Fence_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Fence,%d,%lf\n", ARMCI_Fence_c, ARMCI_Fence_t);
    }

    MPI_Reduce(&ARMCI_Create_mutexes_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Create_mutexes_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Create_mutexes,%d,%lf\n", ARMCI_Create_mutexes_c,
	       ARMCI_Create_mutexes_t);
    }

    MPI_Reduce(&ARMCI_PutS_flag_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutS_flag_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutS_flag,%d,%lf\n", ARMCI_PutS_flag_c,
	       ARMCI_PutS_flag_t);
    }

    MPI_Reduce(&ARMCI_WaitProc_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_WaitProc_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_WaitProc,%d,%lf\n", ARMCI_WaitProc_c,
	       ARMCI_WaitProc_t);
    }

    MPI_Reduce(&ARMCI_Lock_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Lock_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Lock,%d,%lf\n", ARMCI_Lock_c, ARMCI_Lock_t);
    }

    MPI_Reduce(&ARMCI_GetValueDouble_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetValueDouble_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetValueDouble,%d,%lf\n", ARMCI_GetValueDouble_c,
	       ARMCI_GetValueDouble_t);
    }

    MPI_Reduce(&ARMCI_NbGetV_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbGetV_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbGetV,%d,%lf\n", ARMCI_NbGetV_c, ARMCI_NbGetV_t);
    }

    MPI_Reduce(&ARMCI_Rmw_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Rmw_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Rmw,%d,%lf\n", ARMCI_Rmw_c, ARMCI_Rmw_t);
    }

    MPI_Reduce(&ARMCI_Init_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Init_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Init,%d,%lf\n", ARMCI_Init_c, ARMCI_Init_t);
    }

    MPI_Reduce(&ARMCI_WaitAll_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_WaitAll_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_WaitAll,%d,%lf\n", ARMCI_WaitAll_c, ARMCI_WaitAll_t);
    }

    MPI_Reduce(&ARMCI_NbGet_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbGet_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbGet,%d,%lf\n", ARMCI_NbGet_c, ARMCI_NbGet_t);
    }

    MPI_Reduce(&ARMCI_PutValueFloat_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutValueFloat_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutValueFloat,%d,%lf\n", ARMCI_PutValueFloat_c,
	       ARMCI_PutValueFloat_t);
    }

    MPI_Reduce(&ARMCI_NbPutS_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_NbPutS_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_NbPutS,%d,%lf\n", ARMCI_NbPutS_c, ARMCI_NbPutS_t);
    }

    MPI_Reduce(&ARMCI_PutS_flag_dir_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutS_flag_dir_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutS_flag_dir,%d,%lf\n", ARMCI_PutS_flag_dir_c,
	       ARMCI_PutS_flag_dir_t);
    }

    MPI_Reduce(&ARMCI_PutValueLong_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_PutValueLong_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_PutValueLong,%d,%lf\n", ARMCI_PutValueLong_c,
	       ARMCI_PutValueLong_t);
    }

    MPI_Reduce(&ARMCI_Wait_c, &c, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_Wait_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_Wait,%d,%lf\n", ARMCI_Wait_c, ARMCI_Wait_t);
    }

    MPI_Reduce(&ARMCI_GetValueLong_c, &c, 1, MPI_INT, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&ARMCI_GetValueLong_t, &t, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (me == 0) {
	printf("ARMCI_GetValueLong,%d,%lf\n", ARMCI_GetValueLong_c,
	       ARMCI_GetValueLong_t);
    }
    PARMCI_Finalize();
}


int ARMCI_AccV(int op, void *scale, armci_giov_t * darr, int len, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_AccV(op, scale, darr, len, proc);
    etime = TIME();
    ARMCI_AccV_t += etime - stime;
    return rval;
}


void ARMCI_Barrier()
{

    static double stime, etime;
    stime = TIME();
    PARMCI_Barrier();
    etime = TIME();
    ARMCI_Barrier_t += etime - stime;
}


int ARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr,
	       void *dst_ptr, int *dst_stride_arr, int *count,
	       int stride_levels, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_AccS(optype, scale, src_ptr, src_stride_arr, dst_ptr,
		    dst_stride_arr, count, stride_levels, proc);
    etime = TIME();
    ARMCI_AccS_t += etime - stime;
    return rval;
}


int ARMCI_NbPut(void *src, void *dst, int bytes, int proc,
		armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_NbPut(src, dst, bytes, proc, nb_handle);
    etime = TIME();
    ARMCI_NbPut_t += etime - stime;
    return rval;
}


int ARMCI_GetValueInt(void *src, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_GetValueInt(src, proc);
    etime = TIME();
    ARMCI_GetValueInt_t += etime - stime;
    return rval;
}


int ARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v,
		   int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Put_flag(src, dst, bytes, f, v, proc);
    etime = TIME();
    ARMCI_Put_flag_t += etime - stime;
    return rval;
}


int ARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr,
		 int *dst_stride_arr, int *count, int stride_levels,
		 int proc, armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_NbGetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
		      count, stride_levels, proc, nb_handle);
    etime = TIME();
    ARMCI_NbGetS_t += etime - stime;
    return rval;
}


void *ARMCI_Malloc_local(armci_size_t bytes)
{
    void *rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Malloc_local(bytes);
    etime = TIME();
    ARMCI_Malloc_local_t += etime - stime;
    return rval;
}


int ARMCI_Free_local(void *ptr)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Free_local(ptr);
    etime = TIME();
    ARMCI_Free_local_t += etime - stime;
    return rval;
}


int ARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Get(src, dst, bytes, proc);
    etime = TIME();
    ARMCI_Get_t += etime - stime;
    return rval;
}


int ARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Put(src, dst, bytes, proc);
    etime = TIME();
    ARMCI_Put_t += etime - stime;
    return rval;
}


int ARMCI_Destroy_mutexes()
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Destroy_mutexes();
    etime = TIME();
    ARMCI_Destroy_mutexes_t += etime - stime;
    return rval;
}


int ARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr,
	       int *dst_stride_arr, int *count, int stride_levels,
	       int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_GetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
		    count, stride_levels, proc);
    etime = TIME();
    ARMCI_GetS_t += etime - stime;
    return rval;
}


int ARMCI_NbAccV(int op, void *scale, armci_giov_t * darr, int len,
		 int proc, armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_NbAccV(op, scale, darr, len, proc, nb_handle);
    etime = TIME();
    ARMCI_NbAccV_t += etime - stime;
    return rval;
}


float ARMCI_GetValueFloat(void *src, int proc)
{
    float rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_GetValueFloat(src, proc);
    etime = TIME();
    ARMCI_GetValueFloat_t += etime - stime;
    return rval;
}


int ARMCI_Malloc(void **ptr_arr, armci_size_t bytes)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Malloc(ptr_arr, bytes);
    etime = TIME();
    ARMCI_Malloc_t += etime - stime;
    return rval;
}


int ARMCI_NbAccS(int optype, void *scale, void *src_ptr,
		 int *src_stride_arr, void *dst_ptr, int *dst_stride_arr,
		 int *count, int stride_levels, int proc,
		 armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_NbAccS(optype, scale, src_ptr, src_stride_arr, dst_ptr,
		      dst_stride_arr, count, stride_levels, proc,
		      nb_handle);
    etime = TIME();
    ARMCI_NbAccS_t += etime - stime;
    return rval;
}


int ARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr,
	       int *dst_stride_arr, int *count, int stride_levels,
	       int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
		    count, stride_levels, proc);
    etime = TIME();
    ARMCI_PutS_t += etime - stime;
    return rval;
}


int ARMCI_PutV(armci_giov_t * darr, int len, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_PutV(darr, len, proc);
    etime = TIME();
    ARMCI_PutV_t += etime - stime;
    return rval;
}


int ARMCI_Free(void *ptr)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Free(ptr);
    etime = TIME();
    ARMCI_Free_t += etime - stime;
    return rval;
}


int ARMCI_Init_args(int *argc, char ***argv)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Init_args(argc, argv);
    etime = TIME();
    ARMCI_Init_args_t += etime - stime;
    return rval;
}


int ARMCI_PutValueInt(int src, void *dst, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_PutValueInt(src, dst, proc);
    etime = TIME();
    ARMCI_PutValueInt_t += etime - stime;
    return rval;
}


void ARMCI_Memget(size_t bytes, armci_meminfo_t * meminfo, int memflg)
{

    static double stime, etime;
    stime = TIME();
    PARMCI_Memget(bytes, meminfo, memflg);
    etime = TIME();
    ARMCI_Memget_t += etime - stime;
}


void ARMCI_AllFence()
{

    static double stime, etime;
    stime = TIME();
    PARMCI_AllFence();
    etime = TIME();
    ARMCI_AllFence_t += etime - stime;
}


int ARMCI_NbPutV(armci_giov_t * darr, int len, int proc,
		 armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_NbPutV(darr, len, proc, nb_handle);
    etime = TIME();
    ARMCI_NbPutV_t += etime - stime;
    return rval;
}


int ARMCI_PutValueDouble(double src, void *dst, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_PutValueDouble(src, dst, proc);
    etime = TIME();
    ARMCI_PutValueDouble_t += etime - stime;
    return rval;
}


int ARMCI_GetV(armci_giov_t * darr, int len, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_GetV(darr, len, proc);
    etime = TIME();
    ARMCI_GetV_t += etime - stime;
    return rval;
}


int ARMCI_Test(armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Test(nb_handle);
    etime = TIME();
    ARMCI_Test_t += etime - stime;
    return rval;
}


void ARMCI_Unlock(int mutex, int proc)
{

    static double stime, etime;
    stime = TIME();
    PARMCI_Unlock(mutex, proc);
    etime = TIME();
    ARMCI_Unlock_t += etime - stime;
}


void ARMCI_Fence(int proc)
{

    static double stime, etime;
    stime = TIME();
    PARMCI_Fence(proc);
    etime = TIME();
    ARMCI_Fence_t += etime - stime;
}


int ARMCI_Create_mutexes(int num)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Create_mutexes(num);
    etime = TIME();
    ARMCI_Create_mutexes_t += etime - stime;
    return rval;
}


int ARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr,
		    int *dst_stride_arr, int *count, int stride_levels,
		    int *flag, int val, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_PutS_flag(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
			 count, stride_levels, flag, val, proc);
    etime = TIME();
    ARMCI_PutS_flag_t += etime - stime;
    return rval;
}


int ARMCI_WaitProc(int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_WaitProc(proc);
    etime = TIME();
    ARMCI_WaitProc_t += etime - stime;
    return rval;
}


void ARMCI_Lock(int mutex, int proc)
{

    static double stime, etime;
    stime = TIME();
    PARMCI_Lock(mutex, proc);
    etime = TIME();
    ARMCI_Lock_t += etime - stime;
}


double ARMCI_GetValueDouble(void *src, int proc)
{
    double rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_GetValueDouble(src, proc);
    etime = TIME();
    ARMCI_GetValueDouble_t += etime - stime;
    return rval;
}


int ARMCI_NbGetV(armci_giov_t * darr, int len, int proc,
		 armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_NbGetV(darr, len, proc, nb_handle);
    etime = TIME();
    ARMCI_NbGetV_t += etime - stime;
    return rval;
}


int ARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Rmw(op, ploc, prem, extra, proc);
    etime = TIME();
    ARMCI_Rmw_t += etime - stime;
    return rval;
}


int ARMCI_Init()
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Init();
    etime = TIME();
    ARMCI_Init_t += etime - stime;
    return rval;
}


int ARMCI_WaitAll()
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_WaitAll();
    etime = TIME();
    ARMCI_WaitAll_t += etime - stime;
    return rval;
}


int ARMCI_NbGet(void *src, void *dst, int bytes, int proc,
		armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_NbGet(src, dst, bytes, proc, nb_handle);
    etime = TIME();
    ARMCI_NbGet_t += etime - stime;
    return rval;
}


int ARMCI_PutValueFloat(float src, void *dst, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_PutValueFloat(src, dst, proc);
    etime = TIME();
    ARMCI_PutValueFloat_t += etime - stime;
    return rval;
}


int ARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr,
		 int *dst_stride_arr, int *count, int stride_levels,
		 int proc, armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_NbPutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
		      count, stride_levels, proc, nb_handle);
    etime = TIME();
    ARMCI_NbPutS_t += etime - stime;
    return rval;
}


int ARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr,
			int *dst_stride_arr, int *count, int stride_levels,
			int *flag, int val, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval =
	PARMCI_PutS_flag_dir(src_ptr, src_stride_arr, dst_ptr,
			     dst_stride_arr, count, stride_levels, flag,
			     val, proc);
    etime = TIME();
    ARMCI_PutS_flag_dir_t += etime - stime;
    return rval;
}


int ARMCI_PutValueLong(long src, void *dst, int proc)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_PutValueLong(src, dst, proc);
    etime = TIME();
    ARMCI_PutValueLong_t += etime - stime;
    return rval;
}


int ARMCI_Wait(armci_hdl_t * nb_handle)
{
    int rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_Wait(nb_handle);
    etime = TIME();
    ARMCI_Wait_t += etime - stime;
    return rval;
}


long ARMCI_GetValueLong(void *src, int proc)
{
    long rval;
    static double stime, etime;
    stime = TIME();
    rval = PARMCI_GetValueLong(src, proc);
    etime = TIME();
    ARMCI_GetValueLong_t += etime - stime;
    return rval;
}
