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

int comex_put(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    assert(src && dst && bytes);
    COMEXD_put_nbi(src, dst, bytes, proc);
    COMEXD_waitproc(proc);
    return 0;
}

int comex_get(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    COMEXD_get_nbi(src, dst, bytes, proc);
    COMEXD_waitproc(proc);

    return 0;
}

int comex_acc(int datatype, void *scale, void *src_ptr, void *dst_ptr,
                 int bytes, int proc, comex_group_t group)
{
    double *get_buf = (double *)l_state.acc_buf;
    double *_src_buf = (double *)src_ptr;
    double calc_scale = *(double *)scale;
    int m, limit;


    assert(bytes <= l_state.acc_buf_len);
    assert(datatype == COMEX_ACC_DBL);
    assert(get_buf);

    COMEXD_network_lock(proc);
    comex_get(dst_ptr, get_buf, bytes, proc, group);

    for (m=0, limit=bytes/sizeof(double); m<limit; ++m) {
        if (calc_scale == 1.0) {
            get_buf[m] += _src_buf[m];
        }
        else {
            get_buf[m] += calc_scale * _src_buf[m];
        }
    }

    comex_put(get_buf, dst_ptr, bytes, proc, group);
    COMEXD_network_unlock(proc);
    
    return 0;
}


int comex_nbput(void *src, void *dst, int bytes,
        int proc, comex_group_t group, comex_request_t *hdl)
{
    int rc;
    rc = comex_put(src, dst, bytes, proc, group);
    return rc;
}   
    
int comex_nbget(void *src, void *dst, int bytes,
        int proc, comex_group_t group, comex_request_t *hdl)
{
    int rc;
    rc = comex_get(src, dst, bytes, proc, group);
    return 0;
}

