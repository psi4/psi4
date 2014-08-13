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


#ifndef _CBUF_H
#define _CBUF_H

#include <stdio.h>
#include <pthread.h>
#include <infiniband/verbs.h>
#include <malloc.h>
#include <string.h>

#ifdef _IA64_
#define CBUF_FLAG_TYPE uint64_t
#else
#define CBUF_FLAG_TYPE uint32_t
#endif

#if (defined(RDMA_FAST_PATH) || defined(ADAPTIVE_RDMA_FAST_PATH))

#define CBUF_BUFFER_SIZE (viadev_cbuf_total_size -    \
     2*sizeof(CBUF_FLAG_TYPE))

#else

#define CBUF_BUFFER_SIZE  (viadev_cbuf_total_size)

#endif

/*
 * brief justification for cbuf format:
 * descriptor must be aligned (64 bytes).
 * cbuf size must be multiple of this alignment to allow contiguous allocation
 * descriptor and buffer should be contiguous to allow via implementations that
 * optimize contiguous descriptor/data (? how likely ?)
 * need to be able to store send handle in cbuf so that we can mark sends
 * complete when communication completes. don't want to store
 * it in packet header because we don't always need to send over the network.
 * don't want to store at beginning or between desc and buffer (see above) so
 * store at end. 
 */

#define QP_CON_REQ 1
#define QP_CON_ACK 2
#define QP_CON_BREAK_FROM_CLIENT 3
#define QP_CON_BREAK_FROM_SERVER 4

typedef struct {
    uint32_t src_rank;
    uint16_t lid;
    uint32_t qpnum;
    int msg_type;
} armci_ud_rank;

struct ibv_descriptor {
    union {
        struct ibv_recv_wr rr;
        struct ibv_send_wr sr;
    } u;
    struct ibv_sge sg_entry;
    void *next;
};

typedef struct ibv_descriptor IBV_DESCRIPTOR;

typedef struct cbuf {

#if (defined(RDMA_FAST_PATH) || defined(ADAPTIVE_RDMA_FAST_PATH))
    CBUF_FLAG_TYPE *head_flag;
#endif

    unsigned char *buffer;

#if (defined(RDMA_FAST_PATH) || defined(ADAPTIVE_RDMA_FAST_PATH))
    int padding;
#endif

    IBV_DESCRIPTOR desc;

    /* NULL shandle means not send or not complete. Non-null
     * means pointer to send handle that is now complete. Used
     * by viadev_process_send
     */
    void *shandle;
    struct cbuf_region *region;
    int grank;

    uint16_t bytes_remaining;
    uint16_t len;
    unsigned char *data_start;
    uint16_t ref_count;

} cbuf;
/* one for head and one for tail */
#define CBUF_FAST_RDMA_EXTRA_BYTES (2*sizeof(CBUF_FLAG_TYPE))

#define FAST_RDMA_ALT_TAG 0x8000
#define FAST_RDMA_SIZE_MASK 0x7fff

/*
 * Vbufs are allocated in blocks and threaded on a single free list.
 *
 * These data structures record information on all the cbuf
 * regions that have been allocated.  They can be used for
 * error checking and to un-register and deallocate the regions
 * at program termination.
 *
 */
typedef struct cbuf_region {
    struct ibv_mr *mem_handle;  /* memory handle for entire region */

    void *malloc_start;         /* used to free region later */
    void *malloc_end;           /* to bracket mem region */
    void *malloc_buf_start;     /* used to free DMA region later */
    void *malloc_buf_end;       /* bracket DMA region */
    int count;                  /* number of cbufs in region */
    struct cbuf *cbuf_head;     /* first cbuf in region */
    struct cbuf_region *next;   /* thread cbuf regions */
} cbuf_region;

/* ack. needed here after cbuf is defined. need to clean up header files
 * and dependencies.
 */

typedef unsigned long aint_t;

extern int viadev_cbuf_max;
extern int viadev_cbuf_total_size;
extern int viadev_cbuf_secondary_pool_size;


void dump_cbuf_regions(void);
void dump_cbuf_region(cbuf_region * r);


void allocate_cbufs(int ncbufs);
void deallocate_cbufs(void);

cbuf *get_cbuf(void);
void release_cbuf(cbuf * v);
void cbuf_init_send(cbuf * v, unsigned long len);
void cbuf_init_recv(cbuf * v, unsigned long len);
void cbuf_init_sendrecv(cbuf * v, unsigned long len);

void cbuf_init_rput(cbuf * v,
                    void *local_address,
                    uint32_t local_memhandle,
                    void *remote_address,
                    uint32_t remote_memhandle_rkey, int nbytes);
void cbuf_init_rget(cbuf * v,
                    void *local_address,
                    uint32_t local_memhandle,
                    void *remote_address,
                    uint32_t remote_memhandle_rkey, int nbytes);

void init_cbuf_lock();

void dump_cbuf(char *msg, cbuf * v);

struct ibv_mr * armci_register_memory(void *, int);

/* 
 * Macros for working with cbufs 
 */

#define CBUF_BUFFER_START(v) (v->buffer)

#define CBUF_DATA_SIZE(type) (CBUF_BUFFER_SIZE - sizeof(type))

#endif                          /* _CBUF_H */
