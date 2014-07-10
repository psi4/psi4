#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_MALLOC_H
#   include <malloc.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "typesf2c.h"
#include "sf.h"
#include "sff2c.h"

int SF_Create(char* fname, SFsize_t size_hard_limit,
        SFsize_t size_soft_limit, SFsize_t req_size, int *handle)
{
    Integer s_a;
    int retcode = sfi_create(fname, &size_hard_limit, &size_soft_limit,
            &req_size, &s_a);
    *handle = s_a;
    return retcode;
}

 
int SF_Create_suffix(char* fname, SFsize_t size_hard_limit,
        SFsize_t size_soft_limit, SFsize_t req_size, int *handle, int suffix)
{
    Integer s_a;
    Integer isuffix = suffix;
    int retcode = sfi_create_suffix(fname, &size_hard_limit, &size_soft_limit,
            &req_size, &s_a, &isuffix);
    *handle = s_a;
    return retcode;
}


int SF_Destroy(int handle)
{
    Integer s_a = handle;
    return sf_destroy_(&s_a);
}


int SF_Rwtor(int handle)
{
    Integer s_a = handle;
    return sf_rwtor_(&s_a);
}


int SF_Open(int handle)
{
    Integer s_a = handle;
    return sf_open_(&s_a);
}


int SF_Close(int handle)
{
    Integer s_a = handle;
    return sf_close_(&s_a);
}

int SF_Fsync(int handle)
{
    Integer s_a = handle;
    return sf_fsync_(&s_a);
}


int SF_Write(int handle, SFsize_t offset, SFsize_t bytes,
        char *buffer, int *req_id)
{
    Integer s_a = handle;
    Integer id;
    int retcode = sf_write_(&s_a, &offset, &bytes, buffer, &id);
    *req_id = id;
    return retcode;
}


int SF_Read(int handle, SFsize_t offset, SFsize_t bytes,
        char *buffer, int *req_id)
{
    Integer s_a = handle;
    Integer id;
    int retcode = sf_read_(&s_a, &offset, &bytes, buffer, &id);
    *req_id = id;
    return retcode;
}


int SF_Wait(int *req_id)
{
    Integer id = (Integer)*req_id;
    int retcode = sf_wait_(&id);
    *req_id = id;
    return retcode;
}


int SF_Waitall(int *list, int num)
{
    Integer *copy = (Integer*)malloc(num * sizeof(Integer));
    Integer inum = num;
    int retcode;
    int i;
    for (i=0; i<num; ++i) {
        copy[i] = list[i];
    }
    retcode = sf_waitall_(copy, &inum);
    for (i=0; i<num; ++i) {
        list[i] = copy[i];
    }
    free(copy);
    return retcode;
}


void SF_Errmsg(int code, char *msg)
{
    sfi_errmsg(code, msg);
}
