/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/
/********************************************************
 * An example source module to accompany...
 *
 * "Using POSIX Threads: Programming with Pthreads"
 *     by Brad nichols, Dick Buttlar, Jackie Farrell
 *     O'Reilly & Associates, Inc.
 *
 ********************************************************
 * tpool.h --
 *
 * Structures for thread pool
 */

#ifndef _psi_src_bin_detci_tpool_h
#define _psi_src_bin_detci_tpool_h

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

namespace psi { namespace detci {

typedef struct tpool_work {
    void               (*routine)(void *);
    void                *arg;
    struct tpool_work   *next;
} tpool_work_t;

typedef struct tpool {
    /* pool characteristics */
    int                 num_threads;
    int                 max_queue_size;
    int                 do_not_block_when_full;
    /* pool state */
    pthread_t           *threads;
    int                 cur_queue_size;
    tpool_work_t        *queue_head;
    tpool_work_t        *queue_tail;
    int                 queue_closed;
    int                 shutdown;
    int                 threads_awake;
    /* pool synchronization */
    pthread_mutex_t     queue_lock;
    pthread_cond_t      queue_not_empty;
    pthread_cond_t      queue_not_full;
    pthread_cond_t      queue_empty;
    pthread_cond_t      all_work_done;
} *tpool_t;

void tpool_init(
    tpool_t          *tpoolp,
    int              num_threads, 
    int              max_queue_size,
    int              do_not_block_when_full);

int tpool_add_work(
    tpool_t          tpool,
    void             (*routine)(void *),
    void             *arg);

int tpool_destroy(
    tpool_t          tpool,
    int              finish);

// DGAS: This isnt great, cannot run parallel detci as a module
// Will nuke tpool eventually and move to a omp task module
EXTERN tpool_t thread_pool;

void tpool_queue_open(tpool_t tpool);

void tpool_queue_close(tpool_t tpool, int finish);

}} // namespace psi::detci

#endif // header guard