#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <infiniband/verbs.h>
#include "comex.h"
#include "comex_impl.h"
#include "reg_cache.h"

// Registration cache : Defensive Programming

// nprocs: number of processes
// size: number of entries

struct _reg_entry_t **reg_cache;

// cache size for each process
//

int cache_size = -1;

int reg_cache_init(int nprocs, int size)
{

    // Allocate the registration cache:
    reg_cache = (struct _reg_entry_t **)malloc(sizeof(struct _reg_entry_t *) * nprocs); 
    assert(reg_cache); 

    int i;

    for (i = 0; i < nprocs; ++i) {
        reg_cache[i] = (struct _reg_entry_t *)malloc(sizeof(struct _reg_entry_t));
        assert(reg_cache[i]);
        reg_cache[i]->next = NULL;
    }
    return 0;
}


int reg_cache_destroy(int nprocs)
{
    int i;

    // TODO: Deregister all entries
    for (i = 0; i < nprocs; ++i) {
        if (reg_cache[i])
            free(reg_cache[i]);
    }

    assert(reg_cache);
    free(reg_cache);

    return 0;
}

struct _reg_entry_t* reg_cache_find(int rank, void *buf, int len)
{
    struct _reg_entry_t *runner = reg_cache[rank]->next;

    while (runner) {
        if (runner->buf <= buf && runner->len >= len &&
                ((runner->buf + runner->len) >= (buf + len))) {
            return(runner);
        }

        runner = runner->next;
    }

    return NULL;
}

int reg_cache_insert(int rank, void *buf, int len, int lkey, int rkey, struct ibv_mr *mr)
{
    struct _reg_entry_t *node, *runner = (struct _reg_entry_t *)(reg_cache[rank]);

    node = (struct _reg_entry_t *)malloc(sizeof(struct _reg_entry_t));
    assert(node);

    node->buf = buf;
    node->len = len;
    node->lkey = lkey;
    node->rkey = rkey;
    node->next = NULL;

    if (mr) {
        node->mr = mr;
        assert(node->mr);
        assert(node->buf == buf);
    }
    
    
    assert(NULL == reg_cache_find(rank, buf, 0));
    while (runner->next) {
        runner = runner->next;
    }

    runner->next = node;
    return 0;
}

void reg_cache_delete(int rank, void *buf)
{
    struct _reg_entry_t *found, *runner =
        (struct _reg_entry_t *)reg_cache[rank];

    while (runner->next) {
        if (runner->next->buf == buf) {
            found = runner->next;
            runner->next = runner->next->next;
            free(found);
            return;
        }
        runner = runner->next;
    }
    assert(0);

}

#if 0
// Test Driver
int main()
{
    reg_cache_init(16, 4);

    _reg_entry_t *reg_entry;
    long *a = (long *)malloc(sizeof(long) *8192);

    reg_cache_insert(0, a, 16, 1, 1);
    
    reg_entry = reg_cache_find(1, a , 16);
    assert(!reg_entry);
#if 1
    reg_entry = reg_cache_find(0, a , 32);
    assert(!reg_entry);
    
    reg_entry = reg_cache_find(0, a , 16);
    assert(reg_entry);

    reg_entry = reg_cache_find(0, a , 32);
    assert(!reg_entry);
#endif
    reg_cache_destroy(16);
}
#endif
