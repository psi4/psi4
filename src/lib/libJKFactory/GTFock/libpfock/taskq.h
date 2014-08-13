#ifndef __GATASK_H__
#define __GATASK_H__


#include "pfock_def.h"


int init_taskq (PFock_t pfock);

void clean_taskq (PFock_t pfock);

void reset_taskq (PFock_t pfock);

int taskq_next (PFock_t pfock, int myrow, int mycol, int ntasks);


#endif /* __GATASK_H__ */
