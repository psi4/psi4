#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>

#include "comex.h"
#include "comex_impl.h"


int comex_wait_proc(int proc, comex_group_t group)
{   
    return COMEXD_waitproc(proc);
}

int comex_wait(comex_request_t * hdl)
{
    return COMEXD_waitall();
}   

int comex_test(comex_request_t * hdl, int *status)
{
    *status = 0;
    return COMEXD_waitall();
}

int comex_wait_all(comex_group_t group)
{
    return COMEXD_waitall();
}   

