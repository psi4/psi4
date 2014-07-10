/*
 * Copyright (C) 1999-2001 The Regents of the University of California
 * (through E.O. Lawrence Berkeley National Laboratory), subject to
 * approval by the U.S. Department of Energy.
 *
 * Use of this software is under license. The license agreement is included
 * in the file MVICH_LICENSE.TXT.
 *
 * Developed at Berkeley Lab as part of MVICH.
 *
 * Authors: Bill Saphir      <wcsaphir@lbl.gov>
 *          Michael Welcome  <mlwelcome@lbl.gov>
 */

/* Copyright (c) 2002-2008, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT_MVAPICH in the top level MPICH directory.
 *
 */

#define _XOPEN_SOURCE 600

#include "cbuf.h"
#include <assert.h>


/*
 * cbufs
 * 
 * cbufs provide system buffers for VMPICH. They are analogous to mbufs
 * in BSD networking. 
 * The primary motivation for cbufs is that implementing MPI on VIA
 * seems to requiring pre-posting a number of fixed-sized buffers. 
 * These buffers must be registered (pinned). Life is easier if 
 * they are all registered at once so there is only one memory 
 * handle. We manage a fixed-size pool of cbufs that are
 * allocated and pinned when a progam starts up. We manage
 * the free cbuf list as a singly linked list. 
 * 
 *  Two different ways to manage the free list as a singly-linked list. 
 *  1. head and tail pointers. Add to tail, remove from head. 
 *  2. only head pointer, treat as a stack. 
 * 
 *  #1 Eliminates contention between adding to list and removing from list 
 *   Lock-free possible?
 * 
 *  #2 Has slightly less overhead when there is no contention, and is more 
 *  likely to produce a cbuf that is already in cache. 
 * 
 *  Currently anticipate that most access near-term will be single-threaded, 
 *  so go with head only.  (#2)
 */

/* head of list of allocated cbuf regions */
static cbuf_region *cbuf_region_head = NULL;

/* 
 * free_cbuf_head is the head of the free list
 */

static cbuf *free_cbuf_head = NULL;

static int cbuf_n_allocated = 0;
static long num_free_cbuf = 0;
static long num_cbuf_get = 0;
static long num_cbuf_free = 0;

static pthread_spinlock_t cbuf_lock;
int viadev_cbuf_max = -1;
int viadev_cbuf_total_size = (2 * 1024);
int viadev_cbuf_secondary_pool_size = 128;

void init_cbuf_lock()
{
    pthread_spin_init(&cbuf_lock, 0);
}

static void lock_cbuf()
{
    pthread_spin_lock(&cbuf_lock);
    return;
}

static void unlock_cbuf()
{
    pthread_spin_unlock(&cbuf_lock);
    return;
}


void dump_cbuf_region(cbuf_region * r)
{
}

void dump_cbuf_regions()
{
    cbuf_region *r = cbuf_region_head;

    while (r) {
        dump_cbuf_region(r);
        r = r->next;
    }
}
void deallocate_cbufs()
{
    cbuf_region *r = cbuf_region_head;

    lock_cbuf();

    while (r) {
        if (r->mem_handle != NULL) {
            /* free cbufs add it later */
        }
        r = r->next;
    }

    unlock_cbuf();
}

static void allocate_cbuf_region(int ncbufs)
{
    struct cbuf_region *reg;
    void *mem;
    void *cbuf_dma_buffer;

    int i;
    cbuf *cur;
    int alignment_cbuf = 64;
    int alignment_dma;

    alignment_dma = getpagesize();

    if (free_cbuf_head != NULL) {
    }

    if (ncbufs <= 0) {
    }
    
    /* are we limiting cbuf allocation?  If so, make sure
     * we dont alloc more than allowed
     */

    reg = (struct cbuf_region *) malloc(sizeof(struct cbuf_region));
    if (NULL == reg) {
    }

    if(posix_memalign((void **) &mem, alignment_cbuf, ncbufs * sizeof(cbuf))) {
    }

    /* ALLOCATE THE DMA BUFFER */

    if(posix_memalign((void **) &cbuf_dma_buffer, alignment_dma,
                ncbufs * viadev_cbuf_total_size)) {
    }

    memset(mem, 0, ncbufs * sizeof(cbuf));
    memset(cbuf_dma_buffer, 0, ncbufs * viadev_cbuf_total_size);

    cbuf_n_allocated += ncbufs;
    num_free_cbuf += ncbufs;
    reg->malloc_start = mem;

    reg->malloc_buf_start = cbuf_dma_buffer;
    reg->malloc_end = (void *) ((char *) mem + ncbufs * sizeof(cbuf));
    reg->malloc_buf_end = (void *) ((char *) cbuf_dma_buffer + 
            ncbufs * viadev_cbuf_total_size);
    
    reg->count = ncbufs;

    free_cbuf_head = (cbuf *) ((aint_t) mem);
    
    reg->cbuf_head = free_cbuf_head;


    reg->mem_handle = armci_register_memory(cbuf_dma_buffer,
            ncbufs * viadev_cbuf_total_size);

    if (reg->mem_handle == NULL) {
    }

    /* init the free list */
    for (i = 0; i < ncbufs - 1; i++) {
        cur = free_cbuf_head + i;

        cur->desc.next = free_cbuf_head + i + 1;
        cur->region = reg;

#ifdef ADAPTIVE_RDMA_FAST_PATH
#else
        cur->buffer = (unsigned char *) ((char *)(cbuf_dma_buffer) +
                (i * viadev_cbuf_total_size));
#endif
    }
    /* last one needs to be set to NULL */
    cur = free_cbuf_head + ncbufs - 1;

    cur->desc.next = NULL;

    cur->region = reg;

#ifdef ADAPTIVE_RDMA_FAST_PATH
#else
    cur->buffer = (unsigned char *) ((char *)cbuf_dma_buffer + 
            ((ncbufs - 1) * viadev_cbuf_total_size));

#endif

    /* thread region list */
    reg->next = cbuf_region_head;
    cbuf_region_head = reg;

}
void allocate_cbufs(int ncbufs)
{
    /* this function is only called by the init routines.
     * cache the nic handle and ptag for later cbuf_region allocations
     */
    /* now allocate the first cbuf region */
    allocate_cbuf_region(ncbufs);
}


/* 
 * Get a cbuf off the free list 
 */

cbuf *get_cbuf(void)
{
    cbuf *v;

    lock_cbuf();

    /* 
     * It will often be possible for higher layers to recover
     * when no cbuf is available, but waiting for more descriptors
     * to complete. For now, just abort. 
     */
    if (NULL == free_cbuf_head) {
        allocate_cbuf_region(viadev_cbuf_secondary_pool_size);
        if (NULL == free_cbuf_head) {
        }
    }
    v = free_cbuf_head;
    num_free_cbuf--;
    num_cbuf_get++;

    /* this correctly handles removing from single entry free list */
    free_cbuf_head = free_cbuf_head->desc.next;
#ifdef ADAPTIVE_RDMA_FAST_PATH
    /* need to change this to RPUT_CBUF_FLAG or RGET_CBUF_FLAG later
     * if we are doing rput */
    v->padding = NORMAL_CBUF_FLAG;
#endif

    /* this is probably not the right place to initialize shandle to NULL.
     * Do it here for now because it will make sure it is always initialized.
     * Otherwise we would need to very carefully add the initialization in
     * a dozen other places, and probably miss one. 
     */
    v->shandle = NULL;

    v->ref_count = 0;
    v->len = 0;

    v->grank = -1; /* Make sure it is not inadvertantly used anywhere */

    unlock_cbuf();

    return (v);
}

/*
 * Put a cbuf back on the free list 
 */

void release_cbuf(cbuf * v)
{

    lock_cbuf();

    /* note this correctly handles appending to empty free list */


    assert(v != free_cbuf_head);

    v->desc.next = free_cbuf_head;

#ifdef ADAPTIVE_RDMA_FAST_PATH
#endif


    free_cbuf_head = v;
    num_free_cbuf++;
    num_cbuf_free++;

    unlock_cbuf();
}


/*
 * fill in cbuf descriptor with all necessary info 
 */



void cbuf_init_send(cbuf * v, unsigned long len)
{
    v->desc.u.sr.next = NULL;
    v->desc.u.sr.send_flags = IBV_SEND_SIGNALED;
    v->desc.u.sr.opcode = IBV_WR_SEND;
    v->desc.u.sr.wr_id = (aint_t) v;
    v->desc.u.sr.num_sge = 1;
    v->desc.u.sr.sg_list = &(v->desc.sg_entry);

    v->desc.sg_entry.addr = (uintptr_t) v->buffer;
    v->desc.sg_entry.length = len;
    v->desc.sg_entry.lkey = v->region->mem_handle->lkey;
}

void cbuf_init_recv(cbuf * v, unsigned long len)
{
    v->desc.u.rr.next = NULL;
    v->desc.u.rr.wr_id = (aint_t) v;
    v->desc.u.rr.num_sge = 1;
    v->desc.u.rr.sg_list = &(v->desc.sg_entry);

    v->desc.sg_entry.addr = (uintptr_t) v->buffer;
    v->desc.sg_entry.length = len;
    v->desc.sg_entry.lkey = v->region->mem_handle->lkey;

#ifdef ADAPTIVE_RDMA_FAST_PATH
    v->padding = NORMAL_CBUF_FLAG;
#endif
}
void cbuf_init_sendrecv(cbuf * v, unsigned long len)
{
}

void cbuf_init_rput(cbuf * v, void *local_address,
                    uint32_t lkey, void *remote_address,
                    uint32_t rkey, int len)
{
    v->desc.u.sr.next = NULL;
    v->desc.u.sr.send_flags = IBV_SEND_SIGNALED;
    v->desc.u.sr.opcode = IBV_WR_RDMA_WRITE;
    v->desc.u.sr.wr_id = (aint_t) v;

    v->desc.u.sr.num_sge = 1;
    v->desc.u.sr.sg_list = &(v->desc.sg_entry);

    v->desc.sg_entry.length = len;
    v->desc.sg_entry.lkey = lkey;
    v->desc.sg_entry.addr = (uintptr_t) local_address;

    v->desc.u.sr.wr.rdma.remote_addr = (uintptr_t) remote_address;
    v->desc.u.sr.wr.rdma.rkey = rkey;

#ifdef ADAPTIVE_RDMA_FAST_PATH
    v->padding = RPUT_CBUF_FLAG;
#endif

}



void cbuf_init_rget(cbuf * v,
                    void *local_address,
                    uint32_t lkey,
                    void *remote_address,
                    uint32_t rkey, int len)
{
    v->desc.u.sr.next = NULL;
    v->desc.u.sr.send_flags = IBV_SEND_SIGNALED;
    v->desc.u.sr.opcode = IBV_WR_RDMA_READ;
    v->desc.u.sr.wr_id = (aint_t) v;

    v->desc.u.sr.num_sge = 1;
    v->desc.u.sr.sg_list = &(v->desc.sg_entry);

    v->desc.sg_entry.length = len;
    v->desc.sg_entry.lkey = lkey;
    v->desc.sg_entry.addr = (uintptr_t) local_address;

    v->desc.u.sr.wr.rdma.remote_addr = (uintptr_t) remote_address;
    v->desc.u.sr.wr.rdma.rkey = rkey;

#ifdef ADAPTIVE_RDMA_FAST_PATH
    v->padding = RGET_CBUF_FLAG;
#endif

}

/*
 * print out cbuf contents for debugging 
 */

void dump_cbuf(char *msg, cbuf * v)
{
}

#ifdef ADAPTIVE_RDMA_FAST_PATH
#endif
