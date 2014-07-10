#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include <mpi.h>

#include "ga.h"
#include "ga-mpi.h"
#include "macdecls.h"

#define n 4
#define PROC_LIST_SIZE 100

int main(int argc, char **argv) {
    int me;
    int g_a;
    int status;
    int i,j;
    int dims[] = {n,n};
    int proc_group[PROC_LIST_SIZE],proclist[PROC_LIST_SIZE],inode;
    int sbuf[1],rbuf[1];
    MPI_Comm comm;

    MPI_Init(&argc, &argv);
    GA_Initialize();
    me = GA_Nodeid();

    status = MA_init(MT_DBL, 100000, 100000);
    if (!status) GA_Error("ma_init failed",-1);
    status = MA_set_auto_verify(1);
    status = MA_set_hard_fail(1);
    status = MA_set_error_print(1);

    inode = GA_Cluster_nodeid();
    if (me == 0) {
        printf("there are %d nodes, node 0 has %d procs\n",
                GA_Cluster_nnodes(), GA_Cluster_nprocs(0));
        fflush(stdout);
    }
    GA_Sync();
    for (i=0; i<GA_Cluster_nnodes(); ++i) {
        for (j=0; j<GA_Cluster_nprocs(i); ++j) {
            proclist[j]=GA_Cluster_procid(i,j);
        }
        proc_group[i]=GA_Pgroup_create(proclist,GA_Cluster_nprocs(i));
    }
    GA_Sync();
    for (i=0; i<GA_Cluster_nnodes(); ++i) {
        if (i == inode) {
            printf("%d joining group %d\n", me, proc_group[inode]);
            GA_Pgroup_set_default(proc_group[inode]);
            g_a = NGA_Create(C_DBL, 2, dims, "a", NULL);
            if (!g_a) GA_Error("NGA_Create failed",-1);
            printf("%d Created array of  group %d as proc no. %d\n",
                    me, proc_group[inode], GA_Nodeid());
            GA_Print_distribution(g_a);
            comm = GA_MPI_Comm_pgroup_default();
            if (comm != MPI_COMM_NULL) {
                sbuf[0] = GA_Nodeid();
                status = MPI_Allreduce(sbuf, rbuf, 1, MPI_INT, MPI_MAX, comm);
                printf("%d max nodeid is %d\n", me, rbuf[0]);
                if ((rbuf[0]+1) != GA_Cluster_nprocs(i)) {
                    GA_Error("MPI_Allreduce failed",1);
                }
            }
            else {
                printf("MPI_Comm was null!\n");
            }
            GA_Pgroup_set_default(GA_Pgroup_get_world());
        }
        GA_Sync();
    }

    GA_Terminate();
    MPI_Finalize();

    return 0;
}
