#ifndef _REG_CACHE_H_
#define _REG_CACHE_H_

#include <dmapp.h>

/**
 * Enumerate the return codes for registration cache functions.
 */
typedef enum _reg_return_t {
    RR_SUCCESS=0,   /**< success */
    RR_FAILURE      /**< non-specific failure */
} reg_return_t;

/**
 * A registered contiguous memory region.
 */
typedef struct _reg_entry_t {
    void *buf;                  /**< starting address of region */
    size_t len;                 /**< length of region */
    dmapp_seg_desc_t mr;        /**< dmapp registered memory region */
    struct _reg_entry_t *next;  /**< next memory region in list */
} reg_entry_t;

/* functions
 *
 * documentation is in the *.c file
 */

reg_return_t reg_cache_init(int nprocs);
reg_return_t reg_cache_destroy();
reg_entry_t *reg_cache_find(int rank, void *buf, int len);
reg_entry_t *reg_cache_insert(int rank, void *buf, int len, dmapp_seg_desc_t mr);
reg_return_t reg_cache_delete(int, void *buf);

#endif /* _REG_CACHE_H_ */
