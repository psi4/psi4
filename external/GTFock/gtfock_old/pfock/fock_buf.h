#ifndef __FOCK_BUF_H__
#define __FOCK_BUF_H__


#include "pfock.h"
#include "CInt.h"


#define PLEN     3
#define P_LO     0
#define P_HI     1
#define P_W      2


void load_local_bufD(PFock_t pfock);

void store_local_bufF(PFock_t pfock);

void compute_FD_ptr(PFock_t pfock, int startM, int endM,
                    int *ptrrow, int *rowsize);

void init_FD_load(PFock_t pfock, int *ptrrow,
                  int **loadrow, int *loadsize);


#endif /* #define __FOCK_BUF_H__ */
