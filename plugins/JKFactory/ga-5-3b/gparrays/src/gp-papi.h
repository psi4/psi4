#ifndef GPPAPI_H
#define GPPAPI_H

#include "typesf2c.h"
#include <stdlib.h>

/* Routines from gpbase.c */

extern void    pgp_debug(Integer g_p, Integer intsize);
extern void    pgp_initialize();
extern void    pgp_terminate();
extern void*   pgp_malloc(size_t size);
extern void    pgp_free(void* ptr);
extern Integer pgp_create_handle();
extern void    pgp_set_dimensions(Integer g_p, Integer ndim, Integer *dims,
                                  Integer intsize);
extern Integer pgp_get_dimension(Integer g_p);
extern void    pgp_set_chunk(Integer g_p, Integer *chunk);
extern void    pgp_set_irreg_distr(Integer g_p, Integer *mapc, Integer *nblock);
extern logical pgp_allocate(Integer g_p);
extern logical pgp_destroy(Integer g_p);
extern void    pgp_distribution(Integer g_p, Integer proc,
                                Integer *lo, Integer *hi);
extern void    pgp_assign_local_element(Integer g_p, Integer *subscript,
                                        void *ptr, Integer size, Integer intsize);
extern void*   pgp_free_local_element(Integer g_p, Integer *subscript);
extern void    pgp_memzero(Integer g_p, Integer intsize);
extern void    pgp_sync();

/* Routines from gponesided.c */
extern void    pgp_get(Integer g_p, Integer *lo, Integer *hi, void *buf,
                       void **buf_ptr, Integer *ld,
                       void *buf_size, Integer *ld_sz, 
                       Integer *size, Integer isize, Integer setbuf);
extern void    pgp_put(Integer g_p, Integer *lo, Integer *hi, void **buf_ptr,
                       Integer *ld, void *buf_size, Integer *ld_sz, 
                       Integer *size, Integer checksize, Integer isize);
extern void    pgp_get_size(Integer g_p, Integer *lo, Integer *hi,
                            Integer *size, Integer isize);
extern void    pgp_access_element(Integer g_p, Integer *subscript,
                                  void *ptr, Integer *size);
extern void    pgp_release_element(Integer g_p, Integer *subscript);
extern void    pgp_release_update_element(Integer g_p, Integer *subscript);
extern void    pgp_gather_size(Integer g_p, Integer nv, Integer *subscript,
                               Integer *size, Integer intsize);
extern void    pgp_gather(Integer g_p, Integer nv, Integer *subscript,
                          void *buf, void **buf_ptr, void *buf_size,
                          Integer *size, Integer intsize, Integer setbuf);
extern void    pgp_scatter(Integer g_p, Integer nv, Integer *subscript,
                           void **buf_ptr, void *buf_size, Integer *size,
                           Integer checksize, Integer intsize);
#endif  /* GPPAPI_H */
