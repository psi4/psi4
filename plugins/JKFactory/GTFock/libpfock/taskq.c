#include <stdlib.h>
#include <stdio.h>
#include <ga.h>

#include "config.h"
#include "taskq.h"


int init_taskq (PFock_t pfock)
{
    int dims[2];
    int block[2];
    int *map;
    int nprow;
    int npcol;
    int i;

    // create GA for dynamic scheduler
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    dims[0] = nprow;
    dims[1] = npcol;
    block[0] = nprow;
    block[1] = npcol;    
    map = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + npcol));
    if (NULL == map)
    {
        return -1;
    }    
    for (i = 0; i < pfock->nprow; i++)
    {
        map[i] = i;
    }
    for (i = 0; i < npcol; i++)
    {
        map[i + nprow] = i;
    }
    pfock->ga_taskid = NGA_Create_irreg (C_INT, 2, dims, "array taskid", block, map);
    if (0 == pfock->ga_taskid)
    {
        return -1;
    }
    free (map);
    
    return 0;
}


void clean_taskq (PFock_t pfock)
{
    GA_Destroy (pfock->ga_taskid);
}

//RMR removed inline because it was giving problems because the declaration was not in the header file
void reset_taskq (PFock_t pfock)
{
    int izero = 0;    
    GA_Fill (pfock->ga_taskid, &izero);
}

//RMR removed inline because it was giving problems because the declaration was not in the header file
int taskq_next (PFock_t pfock, int myrow, int mycol, int ntasks)
{
    int nxtask;
    int idx[2];

    idx[0] = myrow;
    idx[1] = mycol;
    nxtask = NGA_Read_inc (pfock->ga_taskid, idx, ntasks);   
    return nxtask;
}
