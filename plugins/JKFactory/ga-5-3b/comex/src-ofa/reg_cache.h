#ifndef _REG_CACHE_H_
#define _REG_CACHE_H_

#include <infiniband/verbs.h>

struct _reg_entry_t {
    void *buf;
    int len;
    int lkey;
    int rkey;
    struct ibv_mr *mr;
    struct _reg_entry_t *next;
};

struct _reg_entry_t *reg_cache_find(int, void *, int);
int reg_cache_init(int, int);
int reg_cache_destroy(int);
int reg_cache_insert(int rank, void *buf, int len, int, int, struct ibv_mr *);
void reg_cache_delete(int rank, void *buf);

#endif /* _REG_CACHE_H_ */
