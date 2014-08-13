/**
 * Registration cache.
 *
 * Defensive programming via copious assert statements is encouraged.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C headers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* 3rd party headers */
#include <mpi.h>
#include <dmapp.h>

/* our headers */
#include "comex.h"
#include "comex_impl.h"
#include "reg_cache.h"


/**
 * A registered dmapp segment.
 *
 * dmapp segment registrations *can* return the same segment for instance if
 * the user is allocating many small buffers which are smaller than the page
 * size used by dmapp. We only keep track of local dmapp memory registrations.
 */
typedef struct _dmapp_entry_t {
    dmapp_seg_desc_t mr;            /**< dmapp registered memory region */
    int count;                      /**< ref count */
    struct _dmapp_entry_t *next;    /**< next memory region in list */
} dmapp_entry_t;


/* the static members in this module */
static reg_entry_t **reg_cache = NULL; /**< list of caches (one per process) */
static int reg_nprocs = 0; /**< number of caches (one per process) */
static dmapp_entry_t *dmapp_cache = NULL; /**< list of cached dmapp segments */


/* the static functions in this module */
static dmapp_entry_t *dmapp_cache_find(dmapp_seg_desc_t mr);
static dmapp_entry_t *dmapp_cache_find_intersection(dmapp_seg_desc_t mr);
static dmapp_entry_t *dmapp_cache_insert(dmapp_seg_desc_t mr);
static reg_return_t   dmapp_cache_delete(dmapp_seg_desc_t mr);
static reg_return_t   seg_cmp(void *reg_addr, size_t reg_len,
                              void *oth_addr, size_t oth_len, int op);
static reg_return_t   seg_intersects(void *reg_addr, size_t reg_len,
                                     void *oth_addr, size_t oth_len);
static reg_return_t   seg_contains(void *reg_addr, size_t reg_len,
                                   void *oth_addr, size_t oth_len);
static reg_return_t   reg_entry_intersects(reg_entry_t *reg_entry,
                                           void *buf, int len);
static reg_return_t   reg_entry_contains(reg_entry_t *reg_entry,
                                         void *buf, int len);
static reg_return_t   dmapp_seg_intersects(dmapp_seg_desc_t first,
                                           dmapp_seg_desc_t second);
static reg_return_t   dmapp_seg_contains(dmapp_seg_desc_t first,
                                         dmapp_seg_desc_t second);


#define TEST_FOR_INTERSECTION 0
#define TEST_FOR_CONTAINMENT 1


/**
 * Detects whether two memory segments intersect or one contains the other.
 *
 * @param[in] reg_addr  starting address of original segment
 * @param[in] reg_len   length of original segment
 * @param[in] oth_addr  starting address of other segment
 * @param[in] oth_len   length of other segment
 * @param[in] op        op to perform, either TEST_FOR_INTERSECTION or
 *                      TEST_FOR_CONTAINMENT
 *
 * @pre NULL != reg_beg
 * @pre NULL != oth_beg
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
seg_cmp(void *reg_addr, size_t reg_len, void *oth_addr, size_t oth_len, int op)
{
    ptrdiff_t reg_beg = 0;
    ptrdiff_t reg_end = 0;
    ptrdiff_t oth_beg = 0;
    ptrdiff_t oth_end = 0;
    int result = 0;

    /* preconditions */
    assert(NULL != reg_addr);
    assert(NULL != oth_addr);

    /* casts to ptrdiff_t since arithmetic on void* is undefined */
    reg_beg = (ptrdiff_t)(reg_addr);
    reg_end = reg_beg + (ptrdiff_t)(reg_len) - 1;
    oth_beg = (ptrdiff_t)(oth_addr);
    oth_end = oth_beg + (ptrdiff_t)(oth_len) - 1;

    switch (op) {
        case TEST_FOR_INTERSECTION:
            result = reg_end >= oth_beg && oth_end >= reg_beg;
            break;
        case TEST_FOR_CONTAINMENT:
            result = reg_beg <= oth_beg && oth_end <= reg_end;
            break;
        default:
            assert(0);
    }

    if (result) {
        return RR_SUCCESS;
    }
    else {
        return RR_FAILURE;
    }
}


/**
 * Detects whether two memory segments intersect.
 *
 * @param[in] reg_addr starting address of original segment
 * @param[in] reg_len  length of original segment
 * @param[in] oth_addr starting address of other segment
 * @param[in] oth_len  length of other segment
 *
 * @pre NULL != reg_beg
 * @pre NULL != oth_beg
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
seg_intersects(void *reg_addr, size_t reg_len, void *oth_addr, size_t oth_len)
{
    /* preconditions */
    assert(NULL != reg_addr);
    assert(NULL != oth_addr);

    return seg_cmp(
            reg_addr, reg_len,
            oth_addr, oth_len,
            TEST_FOR_INTERSECTION);
}


/**
 * Detects whether the first memory segment contains the other.
 *
 * @param[in] reg_addr starting address of original segment
 * @param[in] reg_len  length of original segment
 * @param[in] oth_addr starting address of other segment
 * @param[in] oth_len  length of other segment
 *
 * @pre NULL != reg_beg
 * @pre NULL != oth_beg
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
seg_contains(void *reg_addr, size_t reg_len, void *oth_addr, size_t oth_len)
{
    /* preconditions */
    assert(NULL != reg_addr);
    assert(NULL != oth_addr);

    return seg_cmp(
            reg_addr, reg_len,
            oth_addr, oth_len,
            TEST_FOR_CONTAINMENT);
}


/**
 * Detects whether two memory segments intersect.
 *
 * @param[in] reg_entry the registration entry
 * @param[in] buf       starting address for the contiguous memory region
 * @param[in] len       length of the contiguous memory region
 *
 * @pre NULL != reg_entry
 * @pre NULL != buf
 * @pre len >= 0
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
reg_entry_intersects(reg_entry_t *reg_entry, void *buf, int len)
{
    /* preconditions */
    assert(NULL != reg_entry);
    assert(NULL != buf);
    assert(len >= 0);

    return seg_intersects(
            reg_entry->buf, reg_entry->len,
            buf, len);
}


/**
 * Detects whether the first memory segment contains the other.
 *
 * @param[in] reg_entry the registration entry
 * @param[in] buf       starting address for the contiguous memory region
 * @param[in] len       length of the contiguous memory region
 *
 * @pre NULL != reg_entry
 * @pre NULL != buf
 * @pre len >= 0
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
reg_entry_contains(reg_entry_t *reg_entry, void *buf, int len)
{
    /* preconditions */
    assert(NULL != reg_entry);
    assert(NULL != buf);
    assert(len >= 0);

    return seg_contains(
            reg_entry->buf, reg_entry->len,
            buf, len);
}


/**
 * Detects whether two dmapp segments intersect.
 *
 * @param[in] first     the original registration entry
 * @param[in] second    segment to test against
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
dmapp_seg_intersects(dmapp_seg_desc_t first, dmapp_seg_desc_t second)
{
    return seg_intersects(first.addr, first.len, second.addr, second.len);
}


/**
 * Detects whether the first dmapp segment contains the other.
 *
 * @param[in] first     the original registration entry
 * @param[in] second    segment to test against
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
dmapp_seg_contains(dmapp_seg_desc_t first, dmapp_seg_desc_t second)
{
    return seg_contains(first.addr, first.len, second.addr, second.len);
}


/**
 * Remove registration cache entry without deregistration.
 *
 * @param[in] rank the rank where the entry came from
 * @param[in] reg_entry the entry
 *
 * @pre NULL != reg_entry
 * @pre 0 <= rank && rank < reg_nprocs
 *
 * @return RR_SUCCESS on success
 */
static reg_return_t
reg_entry_destroy(int rank, reg_entry_t *reg_entry)
{
    /* preconditions */
    assert(NULL != reg_entry);
    assert(0 <= rank && rank < reg_nprocs);

    if (l_state.rank == rank) {
        dmapp_cache_delete(reg_entry->mr);
    }

    /* free cache entry */
    free(reg_entry);

    return RR_SUCCESS;
}


/**
 * Create internal data structures for the registration cache.
 *
 * @param[in] nprocs    number of registration caches to create i.e. one per
 *                      process
 *
 * @pre this function is called once to initialize the internal data
 * structures and cannot be called again until reg_cache_destroy() has been
 * called
 *
 * @see reg_cache_destroy()
 *
 * @return RR_SUCCESS on success
 */
reg_return_t
reg_cache_init(int nprocs)
{
    int i = 0;

    /* preconditions */
    assert(NULL == reg_cache);
    assert(0 == reg_nprocs);
    assert(NULL == dmapp_cache);

    /* keep the number of caches around for later use */
    reg_nprocs = nprocs;

    /* allocate the registration cache list: */
    reg_cache = (reg_entry_t **)malloc(
            sizeof(reg_entry_t *) * reg_nprocs); 
    assert(reg_cache); 

    /* initialize the registration cache list: */
    for (i = 0; i < reg_nprocs; ++i) {
        reg_cache[i] = NULL;
    }

    return RR_SUCCESS;
}


/**
 * Deregister and destroy all cache entries and associated buffers.
 *
 * @pre this function is called once to destroy the internal data structures
 * and cannot be called again until reg_cache_init() has been called
 *
 * @see reg_cache_init()
 *
 * @return RR_SUCCESS on success
 */
reg_return_t
reg_cache_destroy()
{
    int i = 0;

    /* preconditions */
    assert(NULL != reg_cache);
    assert(0 != reg_nprocs);

    for (i = 0; i < reg_nprocs; ++i) {
        reg_entry_t *runner = reg_cache[i];

        while (runner) {
            reg_entry_t *previous = runner; /* pointer to previous runner */

            /* get next runner */
            runner = runner->next;
            /* destroy the entry */
            reg_entry_destroy(i, previous);
        }
    }

    /* free registration cache list */
    free(reg_cache);
    reg_cache = NULL;

    /* reset the number of caches */
    reg_nprocs = 0;

    /* by the time all entries are destroyed, dmapp cache should be empty */
    assert(NULL == dmapp_cache);

    return RR_SUCCESS;
}


/**
 * Locate a registration cache entry which contains the given segment
 * completely.
 *
 * @param[in] rank  rank of the process
 * @param[in] buf   starting address of the buffer
 * @parma[in] len   length of the buffer
 * 
 * @pre 0 <= rank && rank < reg_nprocs
 * @pre reg_cache_init() was previously called
 *
 * @return the reg cache entry, or NULL on failure
 */
reg_entry_t*
reg_cache_find(int rank, void *buf, int len)
{
    reg_entry_t *entry = NULL;
    reg_entry_t *runner = NULL;

    /* preconditions */
    assert(NULL != reg_cache);
    assert(0 <= rank && rank < reg_nprocs);

    runner = reg_cache[rank];

    while (runner && NULL == entry) {
        if (RR_SUCCESS == reg_entry_contains(runner, buf, len)) {
            entry = runner;
        }
        runner = runner->next;
    }

#ifndef NDEBUG
    /* we assert that the found entry was unique */
    while (runner) {
        if (RR_SUCCESS == reg_entry_contains(runner, buf, len)) {
            assert(0);
        }
        runner = runner->next;
    }
#endif

    return entry;
}


/**
 * Locate a dmapp segment which contains the given segment completely.
 *
 * @param[in] mr    the dmapp segment
 * 
 * @return the reg cache entry, or NULL on failure
 */
static dmapp_entry_t*
dmapp_cache_find(dmapp_seg_desc_t mr)
{
    dmapp_entry_t *entry = NULL;
    dmapp_entry_t *runner = NULL;

    runner = dmapp_cache;

    while (runner && NULL == entry) {
        if (RR_SUCCESS == dmapp_seg_contains(runner->mr, mr)) {
            entry = runner;
        }
        runner = runner->next;
    }

    /* we assert that the found entry was unique */
    while (runner) {
        if (RR_SUCCESS == dmapp_seg_contains(runner->mr, mr)) {
            assert(0);
        }
        runner = runner->next;
    }

    return entry;
}


/**
 * Locate a registration cache entry which intersects the given segment.
 *
 * @param[in] rank  rank of the process
 * @param[in] buf   starting address of the buffer
 * @parma[in] len   length of the buffer
 * 
 * @pre 0 <= rank && rank < reg_nprocs
 * @pre reg_cache_init() was previously called
 *
 * @return the reg cache entry, or NULL on failure
 */
reg_entry_t*
reg_cache_find_intersection(int rank, void *buf, int len)
{
    reg_entry_t *entry = NULL;
    reg_entry_t *runner = NULL;

    /* preconditions */
    assert(NULL != reg_cache);
    assert(0 <= rank && rank < reg_nprocs);

    runner = reg_cache[rank];

    while (runner && NULL == entry) {
        if (RR_SUCCESS == reg_entry_intersects(runner, buf, len)) {
            entry = runner;
        }
        runner = runner->next;
    }

    /* we assert that the found entry was unique */
    while (runner) {
        if (RR_SUCCESS == reg_entry_contains(runner, buf, len)) {
            assert(0);
        }
        runner = runner->next;
    }

    return entry;
}


/**
 * Locate a dmapp segment which intersects the given segment.
 *
 * @param[in] mr    the dmapp segment
 * 
 * @return the reg cache entry, or NULL on failure
 */
static dmapp_entry_t*
dmapp_cache_find_intersection(dmapp_seg_desc_t mr)
{
    dmapp_entry_t *entry = NULL;
    dmapp_entry_t *runner = NULL;

    runner = dmapp_cache;

    while (runner && NULL == entry) {
        if (RR_SUCCESS == dmapp_seg_intersects(runner->mr, mr)) {
            entry = runner;
        }
        runner = runner->next;
    }

    /* we assert that the found entry was unique */
    while (runner) {
        if (RR_SUCCESS == dmapp_seg_intersects(runner->mr, mr)) {
            assert(0);
        }
        runner = runner->next;
    }

    return entry;
}


/**
 * Create a new registration entry based on the given members.
 *
 * @pre 0 <= rank && rank < reg_nprocs
 * @pre NULL != buf
 * @pre 0 <= len
 * @pre reg_cache_init() was previously called
 * @pre NULL == reg_cache_find(rank, buf, 0)
 * @pre NULL == reg_cache_find_intersection(rank, buf, 0)
 *
 * @return RR_SUCCESS on success
 */
reg_entry_t*
reg_cache_insert(int rank, void *buf, int len, dmapp_seg_desc_t mr)
{
    reg_entry_t *node = NULL;

    /* preconditions */
    assert(NULL != reg_cache);
    assert(0 <= rank && rank < reg_nprocs);
    assert(NULL != buf);
    assert(len >= 0);
    assert(NULL == reg_cache_find(rank, buf, len));
    assert(NULL == reg_cache_find_intersection(rank, buf, len));

    if (rank == l_state.rank) {
        dmapp_cache_insert(mr);
    }

    /* allocate the new entry */
    node = (reg_entry_t *)malloc(sizeof(reg_entry_t));
    assert(node);

    /* initialize the new entry */
    node->buf = buf;
    node->len = len;
    node->mr = mr;
    node->next = NULL;

    /* push new entry to tail of linked list */
    if (NULL == reg_cache[rank]) {
        reg_cache[rank] = node;
    }
    else {
        reg_entry_t *runner = reg_cache[rank];
        while (runner->next) {
            runner = runner->next;
        }
        runner->next = node;
    }

    return RR_SUCCESS;
}


/**
 * Removes the reg cache entry associated with the given rank and buffer.
 *
 * If this process owns the buffer, it will unregister the buffer, as well.
 *
 * @param[in] rank
 * @param[in] buf
 *
 * @pre 0 <= rank && rank < reg_nprocs
 * @pre NULL != buf
 * @pre reg_cache_init() was previously called
 * @pre NULL != reg_cache_find(rank, buf, 0)
 *
 * @return RR_SUCCESS on success
 *         RR_FAILURE otherwise
 */
reg_return_t
reg_cache_delete(int rank, void *buf)
{
    reg_return_t status = RR_FAILURE;
    reg_entry_t *runner = NULL;
    reg_entry_t *previous_runner = NULL;

    /* preconditions */
    assert(NULL != reg_cache);
    assert(0 <= rank && rank < reg_nprocs);
    assert(NULL != buf);
    assert(NULL != reg_cache_find(rank, buf, 0));

    /* this is more restrictive than reg_cache_find() in that we locate
     * exactlty the same region starting address */
    runner = reg_cache[rank];
    while (runner) {
        if (runner->buf == buf) {
            break;
        }
        previous_runner = runner;
        runner = runner->next;
    }
    /* we should have found an entry */
    if (NULL == runner) {
        assert(0);
        return RR_FAILURE;
    }

    /* pop the entry out of the linked list */
    if (previous_runner) {
        previous_runner->next = runner->next;
    }
    else {
        reg_cache[rank] = reg_cache[rank]->next;
    }

    status = reg_entry_destroy(rank, runner);

    return status;
}


/**
 * Increments the ref count of an existing dmapp segement or inserts a new
 * entry into the list.
 *
 * @param[in] mr the dmapp segment
 *
 * @return
 */
static dmapp_entry_t *dmapp_cache_insert(dmapp_seg_desc_t mr)
{
    dmapp_entry_t *runner = NULL;
    dmapp_entry_t *previous_runner = NULL;

    /* this is more restrictive than dmapp_cache_find() in that we locate
     * exactlty the same region starting address */
    runner = dmapp_cache;
    while (runner) {
        if (runner->mr.addr == mr.addr) {
            break;
        }
        previous_runner = runner;
        runner = runner->next;
    }

    if (runner) {
        /* make sure it's an exact match */
        assert(runner->mr.len == mr.len);
        /* increment ref count */
        ++(runner->count);
#if DEBUG
        printf("[%d] incrementing ref count of (%p,%zu) to %d\n",
                l_state.rank, runner->mr.addr, runner->mr.len, runner->count);
#endif
    }
    else {
        runner = malloc(sizeof(dmapp_entry_t));
        runner->mr = mr;
        runner->count = 1;
        runner->next = NULL;
        if (previous_runner) {
            previous_runner->next = runner;
        }
        else {
            dmapp_cache = runner;
        }
#if DEBUG
        printf("[%d] inserting (%p,%zu)\n",
                l_state.rank, runner->mr.addr, runner->mr.len);
#endif
    }

    return runner;
}


/**
 * Decrements the ref count and possibly removes the dmapp cache entry
 * associated with the given dmapp segment.
 *
 * @param[in] mr the dmapp segment to possibly remove
 *
 * @pre NULL != dmapp_cache_find(mr)
 *
 * @return RR_SUCCESS on success
 *         RR_FAILURE otherwise
 */
reg_return_t
dmapp_cache_delete(dmapp_seg_desc_t mr)
{
    reg_return_t status = RR_FAILURE;
    dmapp_entry_t *runner = NULL;
    dmapp_entry_t *previous_runner = NULL;

    /* preconditions */
    assert(NULL != dmapp_cache);
    if (NULL == dmapp_cache_find(mr)) {
        printf("[%d] dmapp_cache_find(mr) failed, mr=(%p,%zu)\n",
                l_state.rank, mr.addr, mr.len);
        assert(0);
    }

    /* this is more restrictive than dmapp_cache_find() in that we locate
     * exactlty the same region starting address */
    runner = dmapp_cache;
    while (runner) {
        if (runner->mr.addr == mr.addr) {
            break;
        }
        previous_runner = runner;
        runner = runner->next;
    }
    /* we should have found an entry */
    if (NULL == runner) {
        assert(0);
        return RR_FAILURE;
    }

    /* decrement ref count */
    --(runner->count);
    assert(runner->count >= 0);
#if DEBUG
    printf("[%d] decrementing ref count of (%p,%zu) to %d\n",
            l_state.rank, runner->mr.addr, runner->mr.len, runner->count);
#endif

    if (0 == runner->count) {
        /* pop the entry out of the linked list */
        if (previous_runner) {
            previous_runner->next = runner->next;
        }
        else {
            dmapp_cache = dmapp_cache->next;
        }
#if DEBUG
        printf("[%d] removing (%p,%zu)\n",
                l_state.rank, runner->mr.addr, runner->mr.len);
#endif
        dmapp_mem_unregister(&(runner->mr));

        free(runner);
    }

    return status;
}

