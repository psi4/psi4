#ifndef __FOCKTASK_H__
#define __FOCKTASK_H__


#include "pfock_def.h"
#include "CInt.h"


#define PLEN     3
#define P_LO     0
#define P_HI     1
#define P_W      2


void load_local_bufD (PFock_t pfock);

void store_local_bufF (PFock_t pfock);

void compute_FD_ptr (PFock_t pfock, int startM, int endM,
                     int *ptrrow, int *rowsize);

void init_FD_load (PFock_t pfock, int *ptrrow,
                   int **loadrow, int *loadsize);

void fock_task (PFock_t pfock, BasisSet_t basis, int startrow, int startcol,
                int startM, int endM, int startP, int endP,
                double **D1, double **D2, double **D3,
                double ***VJ1, double ***VJ2, double ***VK3,
                int ldX1, int ldX2, int ldX3,
                double *nsq, double *nitl,double *inttime);

void correct_F (PFock_t pfock);

void access_bufD_GArrays (PFock_t pfock);

void release_bufD_GArrays (PFock_t pfock);


#endif /* #define __FOCKTASK_H__ */
