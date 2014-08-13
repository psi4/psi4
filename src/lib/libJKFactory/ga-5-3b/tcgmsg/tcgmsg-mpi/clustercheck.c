#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_WINSOCK_H
#   include <winsock.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#define HOSTNAME_LEN 128 
static char myname[HOSTNAME_LEN], rootname[HOSTNAME_LEN];

/**
 * return 1 if all processes are running on the same machine, 0 otherwise
 */
int single_cluster()
{
    int rc,me,root=0,stat,global_stat,len;

    gethostname(myname, HOSTNAME_LEN-1);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if(me==root) {
        rc = MPI_Bcast(myname,HOSTNAME_LEN,MPI_CHAR,root,MPI_COMM_WORLD);
    } else {
        rc = MPI_Bcast(rootname, HOSTNAME_LEN,MPI_CHAR,root,MPI_COMM_WORLD);
    }

    if(rc != MPI_SUCCESS){
        fprintf(stderr,"single_cluster:MPI_Bcast failed rc=%d\n",rc);
        MPI_Abort(MPI_COMM_WORLD,rc);
    } 

    len = strlen(myname);
    stat = (me==root) ? 0 : strncmp(rootname, myname, len);

    if(stat != 0) {
        stat = 1;
    }

    rc = MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(rc != MPI_SUCCESS){
        fprintf(stderr,"single_cluster:MPI_MPI_Allreduce failed rc=%d\n",rc);
        MPI_Abort(MPI_COMM_WORLD,rc);
    } 

    if(global_stat) {
        return 0;
    } else {
        return 1;
    }
}
