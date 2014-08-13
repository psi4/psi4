
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "gp-papi.h"
#include "typesf2c.h"


logical wgp_allocate(Integer g_p)
{
    return pgp_allocate(g_p);
}


void wgp_assign_local_element(Integer g_p, Integer *subscript, void *ptr, Integer size)
{
    pgp_assign_local_element(g_p, subscript, ptr, size);
}


Integer wgp_create_handle()
{
    return pgp_create_handle();
}


void wgp_debug(Integer g_p)
{
    pgp_debug(g_p);
}


logical wgp_destroy(Integer g_p)
{
    return pgp_destroy(g_p);
}


void wgp_distribution(Integer g_p, Integer proc, Integer *lo, Integer *hi)
{
    pgp_distribution(g_p, proc, lo, hi);
}


void wgp_free(void *ptr)
{
    pgp_free(ptr);
}


void* wgp_free_local_element(Integer g_p, Integer *subscript)
{
    pgp_free_local_element(g_p, subscript);
}


void wgp_get(Integer g_p, Integer *lo, Integer *hi, void *buf, void **buf_ptr, Integer *ld, void *buf_size, Integer *ld_sz, Integer *size, Integer isize)
{
    pgp_get(g_p, lo, hi, buf, buf_ptr, ld, buf_size, ld_sz, size, isize);
}


Integer wgp_get_dimension(Integer g_p)
{
    return pgp_get_dimension(g_p);
}


void wgp_get_size(Integer g_p, Integer *lo, Integer *hi, Integer *size, Integer isize)
{
    pgp_get_size(g_p, lo, hi, size, isize);
}


void wgp_initialize()
{
    pgp_initialize();
}


void* wgp_malloc(size_t size)
{
    pgp_malloc(size);
}


void wgp_set_chunk(Integer g_p, Integer *chunk)
{
    pgp_set_chunk(g_p, chunk);
}


void wgp_set_dimensions(Integer g_p, Integer ndim, Integer *dims)
{
    pgp_set_dimensions(g_p, ndim, dims);
}


void wgp_terminate()
{
    pgp_terminate();
}

