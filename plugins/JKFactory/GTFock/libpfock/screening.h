#ifndef __SCREENING_H__
#define __SCREENING_H__


#include "pfock_def.h"
#include "CInt.h"


int init_screening (PFock_t pfock, BasisSet_t basis);

int schwartz_screening (PFock_t pfock, int startshell, int endshell);

void clean_screening (PFock_t pfock);


#endif /* __SCREENING_H__ */
