#include <stdlib.h>
#include <stdio.h>
#include <ga.h>

#include "config.h"
#include "taskq.h"


int init_taskq(PFock_t pfock)
{
    int dims[2];
    int block[2];
    
    // create GA for dynamic scheduler
    int nprow = pfock->nprow;
    int npcol = pfock->npcol;
    dims[0] = nprow;
    dims[1] = npcol;
    block[0] = nprow;
    block[1] = npcol;    
    int *map = (int *)PFOCK_MALLOC(sizeof(int) * (nprow + npcol));
    if (NULL == map) {
        return -1;
    }    
    for (int i = 0; i < pfock->nprow; i++) {
        map[i] = i;
    }
    for (int i = 0; i < npcol; i++) {
        map[i + nprow] = i;
    }
    pfock->ga_taskid =
        NGA_Create_irreg(C_INT, 2, dims, "array taskid", block, map);
    if (0 == pfock->ga_taskid) {
        return -1;
    }
    PFOCK_FREE(map);
    
    return 0;
}


void clean_taskq(PFock_t pfock)
{
    GA_Destroy(pfock->ga_taskid);
}


void reset_taskq(PFock_t pfock)
{
    int izero = 0;    
    GA_Fill(pfock->ga_taskid, &izero);
}


int taskq_next(PFock_t pfock, int myrow, int mycol, int ntasks)
{
    int idx[2];

    idx[0] = myrow;
    idx[1] = mycol;
    int nxtask = NGA_Read_inc(pfock->ga_taskid, idx, ntasks);   
    return nxtask;
}