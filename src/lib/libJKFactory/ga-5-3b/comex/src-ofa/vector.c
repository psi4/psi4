#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>

#include "comex.h"
#include "comex_impl.h"

int comex_putv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group)
{
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            assert(src[j] && dst[j] && bytes && limit);
            comex_put(src[j], dst[j], bytes, proc, group);
        }
    }
    return COMEX_SUCCESS;
}

int comex_getv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group)
{
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            comex_get(src[j], dst[j], bytes, proc, group);
        }
    } 
    return COMEX_SUCCESS;
}

int comex_accv(int datatype, void *scale, comex_giov_t *iov,
        int iov_len, int proc, comex_group_t group)
{   
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes; 
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) { 
            comex_acc(datatype, scale, src[j], dst[j], bytes, proc, group);
        }
    }
    return COMEX_SUCCESS;
}   
    
int comex_nbputv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{   
    return comex_putv(iov, iov_len, proc, group);
}

int comex_nbgetv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{   
    return comex_getv(iov, iov_len, proc, group);
}

int comex_nbaccv(int datatype, void *scale, comex_giov_t *iov,
        int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{
    return comex_accv(datatype, scale, iov, iov_len, proc, group);
}
