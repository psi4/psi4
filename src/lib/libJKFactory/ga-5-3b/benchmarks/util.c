#include <stdio.h>
#include <stdlib.h>

#include "ga.h"
#include "typesf2c.h"

#define F77_SET_MA_USE_ARMCI_MEM(name)                                      \
void name()                                                                 \
{                                                                           \
    int retval;                                                             \
    if((retval=setenv("MA_USE_ARMCI_MEM", "YES", 1)) != 0)                  \
        GA_Error("setenv failed: insufficient space in the environment",1); \
}
F77_SET_MA_USE_ARMCI_MEM(set_ma_use_armci_mem)
F77_SET_MA_USE_ARMCI_MEM(set_ma_use_armci_mem_)
F77_SET_MA_USE_ARMCI_MEM(set_ma_use_armci_mem__)
F77_SET_MA_USE_ARMCI_MEM(SET_MA_USE_ARMCI_MEM)
F77_SET_MA_USE_ARMCI_MEM(SET_MA_USE_ARMCI_MEM_)
F77_SET_MA_USE_ARMCI_MEM(SET_MA_USE_ARMCI_MEM__)

#define F77_UTIL_MDTOB(name)                           \
Integer name(Integer *n)                               \
{                                                      \
    if (*n < 0)                                        \
        GA_Error("util_MDTOB_: negative argument",*n); \
    return (Integer) (*n * sizeof(DoublePrecision));   \
}
F77_UTIL_MDTOB(util_mdtob)
F77_UTIL_MDTOB(util_mdtob_)
F77_UTIL_MDTOB(util_mdtob__)
F77_UTIL_MDTOB(UTIL_MDTOB)
F77_UTIL_MDTOB(UTIL_MDTOB_)
F77_UTIL_MDTOB(UTIL_MDTOB__)

#define F77_UTIL_TIMER(name) \
double name()                \
{                            \
    return GA_Wtime();       \
}
F77_UTIL_TIMER(util_timer)
F77_UTIL_TIMER(util_timer_)
F77_UTIL_TIMER(util_timer__)
F77_UTIL_TIMER(UTIL_TIMER)
F77_UTIL_TIMER(UTIL_TIMER_)
F77_UTIL_TIMER(UTIL_TIMER__)

#define F77_FLUSH(name)  \
void name(Integer *unit) \
{                        \
    fflush(stdout);      \
    fflush(stderr);      \
}
F77_FLUSH(ffflush)
F77_FLUSH(ffflush_)
F77_FLUSH(ffflush__)
F77_FLUSH(FFFLUSH)
F77_FLUSH(FFFLUSH_)
F77_FLUSH(FFFLUSH__)
