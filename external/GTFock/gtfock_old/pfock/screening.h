#ifndef __SCREENING_H__
#define __SCREENING_H__


#include "pfock.h"
#include "CInt.h"


int schwartz_screening(PFock_t pfock, BasisSet_t basis);

void clean_screening(PFock_t pfock);


#endif /* __SCREENING_H__ */
