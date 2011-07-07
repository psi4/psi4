/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

#include <world/worldrmi.h>
#include <world/posixmem.h>
#include <world/worldtime.h>
#include <iostream>
#include <algorithm>
#include <utility>

namespace madness {

    RMI* RMI::instance_ptr = 0;

    bool RMI::is_ordered(attrT attr) {
        return attr & ATTR_ORDERED;
    }

    void RMI::run() {
        ThreadBase::set_affinity(1); // The RMI thread is logical thread 1
        if (debugging)
            std::cerr << rank << ":RMI: server thread is running" << std::endl;
        // The RMI server thread spends its life in here
        MPI::Status status[NRECV+1];
        int ind[NRECV+1];
        qmsg q[MAXQ];
        int n_in_q = 0;

        while (1) {

            if (debugging && n_in_q)
                std::cerr << rank << ":RMI: about to call Waitsome with "
                          << n_in_q << " messages in the queue" << std::endl;

            // If MPI is not safe for simultaneous entry by multiple threads we
            // cannot call Waitsome ... have to poll via Testsome
            int narrived;

            MutexWaiter waiter;
            while (!(narrived = SafeMPI::Request::Testsome(NRECV+1, recv_req, ind, status))) {
                if (finished) return;
#if defined(HAVE_CRAYXT) || defined(HAVE_IBMBGP)
                myusleep(1);
#else
                waiter.wait();
#endif
            }

#ifndef HAVE_CRAYXT
            waiter.reset();
#endif

            if (debugging)
                std::cerr << rank << ":RMI: " << narrived
                          << " messages just arrived" << std::endl;

            if (narrived) {
                for (int m=0; m<narrived; ++m) {
                    int src = status[m].Get_source();
                    size_t len = status[m].Get_count(MPI::BYTE);
                    int i = ind[m];

                    ++(stats.nmsg_recv);
                    stats.nbyte_recv += len;

                    const header* h = (const header*)(recv_buf[i]);
                    rmi_handlerT func = h->func;
                    attrT attr = h->attr;
                    counterT count = (attr>>16); //&&0xffff;

                    if (!is_ordered(attr) || count==recv_counters[src]) {
                        // Unordered and in order messages should be digested as soon as possible.
                        if (debugging)
                            std::cerr << rank
                                      << ":RMI: invoking from=" << src
                                      << " nbyte=" << len
                                      << " func=" << func
                                      << " ordered=" << is_ordered(attr)
                                      << " count=" << count
                                      << std::endl;

                        if (is_ordered(attr)) ++(recv_counters[src]);
                        func(recv_buf[i], len);
                        post_recv_buf(i);
                    }
                    else {
                        if (debugging)
                            std::cerr << rank
                                      << ":RMI: enqueing from=" << src
                                      << " nbyte=" << len
                                      << " func=" << func
                                      << " ordered=" << is_ordered(attr)
                                      << " fromcount=" << count
                                      << " herecount=" << int(recv_counters[src])
                                      << std::endl;
                        // Shove it in the queue
                        int n = n_in_q++;
                        if (n >= MAXQ) MADNESS_EXCEPTION("RMI:server: overflowed out-of-order message q\n", n);
                        q[n] = qmsg(len, func, i, src, attr, count);
                    }
                }

                // Only ordered messages can end up in the queue due to
                // out-of-order receipt or order of recv buffer processing.

                // Sort queued messages by ascending recv count
                std::sort(&q[0],&q[0]+n_in_q);

                // Loop thru messages ... since we have sorted only one pass
                // is necessary and if we cannot process a message we
                // save it at the beginning of the queue
                int nleftover = 0;
                for (int m=0; m<n_in_q; ++m) {
                    int src = q[m].src;
                    if (q[m].count == recv_counters[src]) {
                        if (debugging)
                            std::cerr << rank
                                      << ":RMI: queue invoking from=" << src
                                      << " nbyte=" << q[m].len
                                      << " func=" << q[m].func
                                      << " ordered=" << is_ordered(q[m].attr)
                                      << " count=" << q[m].count
                                      << std::endl;

                        ++(recv_counters[src]);
                        q[m].func(recv_buf[q[m].i], q[m].len);
                        post_recv_buf(q[m].i);
                    }
                    else {
                        q[nleftover++] = q[m];
                        if (debugging)
                            std::cerr << rank
                                      << ":RMI: queue pending out of order from=" << src
                                      << " nbyte=" << q[m].len
                                      << " func=" << q[m].func
                                      << " ordered=" << is_ordered(q[m].attr)
                                      << " count=" << q[m].count
                                      << std::endl;
                    }
                }
                n_in_q = nleftover;

                post_pending_huge_msg();
            }
        }
    }

    void RMI::post_pending_huge_msg() {
        if (recv_buf[NRECV]) return;      // Message already pending
        if (!hugeq.empty()) {
            int src = hugeq.front().first;
            size_t nbyte = hugeq.front().second;
            hugeq.pop_front();
            if (posix_memalign((void **)(recv_buf+NRECV), ALIGNMENT, nbyte))
                MADNESS_EXCEPTION("RMI: failed allocating huge message", 1);
            recv_req[NRECV] = comm.Irecv(recv_buf[NRECV], nbyte, MPI::BYTE, src, SafeMPI::RMI_HUGE_DAT_TAG);
            int nada=0;
            comm.Send(&nada, sizeof(nada), MPI::BYTE, src, SafeMPI::RMI_HUGE_ACK_TAG);
        }
    }

    void RMI::post_recv_buf(int i) {
        if (i < NRECV) {
            recv_req[i] = comm.Irecv(recv_buf[i], MAX_MSG_LEN, MPI::BYTE, MPI::ANY_SOURCE, SafeMPI::RMI_TAG);
        }
        else if (i == NRECV) {
            free(recv_buf[i]);
            recv_buf[i] = 0;
            post_pending_huge_msg();
        }
        else {
            MADNESS_EXCEPTION("RMI::post_recv_buf: confusion", i);
        }
    }

    RMI::~RMI() {
        //delete send_counters;
        //delete recv_counters;
        //         if (!MPI::Is_finalized()) {
        //             for (int i=0; i<NRECV; ++i) {
        //                 if (!recv_req[i].Test())
        //                     recv_req[i].Cancel();
        //             }
        //         }
        //for (int i=0; i<NRECV; ++i) free(recv_buf[i]);
    }

    RMI::RMI()
            : comm(MPI::COMM_WORLD)
            , nproc(comm.Get_size())
            , rank(comm.Get_rank())
            , debugging(false)
            , finished(false)
            , send_counters(new unsigned short[nproc])
            , recv_counters(new unsigned short[nproc]) {
        for (int i=0; i<nproc; ++i) send_counters[i] = 0;
        for (int i=0; i<nproc; ++i) recv_counters[i] = 0;
        if (nproc > 1) {
            for (int i=0; i<NRECV; ++i) {
                if (posix_memalign((void**)(recv_buf+i), ALIGNMENT, MAX_MSG_LEN))
                    MADNESS_EXCEPTION("RMI:initialize:failed allocating aligned recv buffer", 1);
                post_recv_buf(i);
            }
            recv_buf[NRECV] = 0;
            start();
        }
    }

    RMI* RMI::instance() {
        if (!instance_ptr) {
            instance_ptr = new RMI();
        }
        return instance_ptr;
    }

    void RMI::huge_msg_handler(void *buf, size_t /*nbytein*/) {
        const size_t* info = (size_t *)(buf);
        int nword = HEADER_LEN/sizeof(size_t);
        int src = info[nword];
        size_t nbyte = info[nword+1];

        instance()->hugeq.push_back(std::make_pair(src,nbyte));
        instance()->post_pending_huge_msg();
    }

    RMI::Request RMI::private_isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, attrT attr) {
        int tag = SafeMPI::RMI_TAG;

        if (nbyte > MAX_MSG_LEN) {
            // Huge message protocol ... send message to dest indicating size and origin of huge message.
            // Remote end posts a buffer then acks the request.  This end can then send.
            const int nword = HEADER_LEN/sizeof(size_t);
            size_t info[nword+2];
            info[nword  ] = rank;
            info[nword+1] = nbyte;

            int ack;
            Request req_ack = comm.Irecv(&ack, sizeof(ack), MPI::BYTE, dest, SafeMPI::RMI_HUGE_ACK_TAG);
            Request req_send = private_isend(info, sizeof(info), dest, RMI::huge_msg_handler, ATTR_UNORDERED);

            MutexWaiter waiter;
            while (!req_send.Test()) waiter.wait();
            waiter.reset();
            while (!req_ack.Test()) waiter.wait();

            tag = SafeMPI::RMI_HUGE_DAT_TAG;
        }
        else if (nbyte < HEADER_LEN) {
            MADNESS_EXCEPTION("RMI::isend --- your buffer is too small to hold the header", static_cast<int>(nbyte));
        }

        if (debugging)
            std::cerr << instance_ptr->rank
                      << ":RMI: sending buf=" << buf
                      << " nbyte=" << nbyte
                      << " dest=" << dest
                      << " func=" << func
                      << " ordered=" << is_ordered(attr)
                      << " count=" << int(send_counters[dest])
                      << std::endl;

        // Since most uses are ordered and we need the mutex to accumulate stats
        // we presently always get the lock
        lock();

        // If ordering need the mutex to enclose sending the message
        // otherwise there is a livelock scenario due to a starved thread
        // holding an early counter.
        if (is_ordered(attr)) {
            //lock();
            attr |= ((send_counters[dest]++)<<16);
        }

        header* h = (header*)(buf);
        h->func = func;
        h->attr = attr;

        ++(stats.nmsg_sent);
        stats.nbyte_sent += nbyte;

        Request result = comm.Isend(buf, nbyte, MPI::BYTE, dest, tag);

        //if (is_ordered(attr)) unlock();
        unlock();

        return result;
    }

    void RMI::private_exit() {
        if (debugging)
            std::cerr << instance_ptr->rank << ":RMI: sending exit request to server thread" << std::endl;

        finished = true;
        myusleep(10000);

        //delete this;
    }

    RMI::Request RMI::isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, unsigned int attr) {
        return instance()->private_isend(buf, nbyte, dest, func, attr);
    }

    void RMI::end() {
        if (instance_ptr) instance_ptr->private_exit();
    }

    void RMI::begin() {
        instance();
    }

    void RMI::set_debug(bool status) {
        instance()->debugging = status;
    }

    bool RMI::get_debug() {
        return instance()->debugging;
    }

    const RMIStats& RMI::get_stats() {
        return instance()->stats;
    }

} // namespace madness
