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

  $Id: worldrmi.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_WORLD_WORLDRMI_H__INCLUDED
#define MADNESS_WORLD_WORLDRMI_H__INCLUDED

#include <world/safempi.h>
//#include <cstdlib>
//#include <iostream>
//#include <algorithm>
#include <world/worldthread.h>
#include <world/worldtypes.h>
//#include <world/posixmem.h>
#include <utility>
#include <list>

/*
  There is just one server thread and it is the only one
  messing with the recv buffers, so there is no need for
  mutex on recv related data.

  Multiple threads (including the server) may send hence
  we need to be careful about send-related data.

  When MPI is initialized we need to use init_thread with
  multiple required.

  This RMI service operates only in COMM_WORLD.  It easy enough
  to extend to other communicators but the point is to have
  only one server thread for all possible uses.  You just
  have to translate rank_in_comm into rank_in_world by
  getting the groups from both communicators using
  MPI_Comm_group and then creating a map from ranks in
  comm to ranks in world using MPI_Group_translate_ranks.

  The class is a singleton ... i.e., there is only one instance of it
  that is made the first time that you call RMI::instance().

  Handler routines should have this type

  typedef void (*rmi_handlerT)(void* buf, size_t nbyte);

  There are few user accessible routines.

  RMI::Request RMI::isend(const void* buf, size_t nbyte, int dest,
                          rmi_handlerT func, unsigned int attr=0)
  - to send an asynchronous message
  - RMI::Request has the same interface as MPI::Request
  (right now it is an MPI::Request but this is not guaranteed)

  void RMI::begin()
  - to start the server thread

  void RMI::end()
  - to terminate the server thread

  bool RMI::get_debug()
  - to get the debug flag

  void RMI::set_debug(bool)
  - to set the debug flag

*/

namespace madness {

    // This is the generic low-level interface for a message handler
    typedef void (*rmi_handlerT)(void* buf, size_t nbyte);

    struct qmsg {
        typedef uint16_t counterT;
        typedef uint32_t attrT;
        size_t len;
        rmi_handlerT func;
        int i;               // buffer index
        ProcessID src;
        attrT attr;
        counterT count;

        qmsg(size_t len, rmi_handlerT func, int i, int src, attrT attr, counterT count)
            : len(len), func(func), i(i), src(src), attr(attr), count(count) {}

        bool operator<(const qmsg& other) const {
            return count < other.count;
        }

        qmsg() {}
    };


    // Holds message passing statistics
    struct RMIStats {
        uint64_t nmsg_sent;
        uint64_t nbyte_sent;
        uint64_t nmsg_recv;
        uint64_t nbyte_recv;

        RMIStats()
                : nmsg_sent(0), nbyte_sent(0), nmsg_recv(0), nbyte_recv(0) {}
    };

    class RMI : private madness::ThreadBase , private madness::Mutex {
        typedef uint16_t counterT;
        typedef uint32_t attrT;
    public:
        // Choose header length to hold at least sizeof(header) and
        // also to ensure good alignment of the user payload.
        static const size_t ALIGNMENT = 64;
        static const size_t HEADER_LEN = ALIGNMENT;
        //static const size_t MAX_MSG_LEN = 256*1024;
        static const size_t MAX_MSG_LEN = 3*512*1024;

        static const attrT ATTR_UNORDERED=0x0;
        static const attrT ATTR_ORDERED=0x1;

        typedef SafeMPI::Request Request;

    private:
#ifdef HAVE_CRAYXT
        static const int NRECV=128;
#else
        static const int NRECV=32;
#endif
        static const int MAXQ=NRECV+1;

        std::list< std::pair<int,size_t> > hugeq; // q for incoming huge messages

        RMIStats stats;
        SafeMPI::Intracomm comm;
        const int nproc;            // No. of processes in comm world
        const ProcessID rank;       // Rank of this process
        volatile bool debugging;    // True if debugging
        volatile bool finished;     // True if finished

        volatile counterT* send_counters;
        counterT* recv_counters;
        unsigned char* recv_buf[NRECV+1]; // Will be at least ALIGNMENT aligned ... +1 for huge messages
        SafeMPI::Request recv_req[NRECV+1];

        static RMI* instance_ptr;    // Pointer to the singleton instance

        static bool is_ordered(attrT attr);

        struct header {
            rmi_handlerT func;
            attrT attr;
        };

        void run();

        void post_pending_huge_msg();

        void post_recv_buf(int i);

        virtual ~RMI();

        RMI();

        static RMI* instance();

        static void huge_msg_handler(void *buf, size_t nbytein);

        Request private_isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, attrT attr);

        void private_exit();

    public:
        static Request isend(const void* buf, size_t nbyte, ProcessID dest, rmi_handlerT func, unsigned int attr=ATTR_UNORDERED);

        static void end();

        static void begin();

        static void set_debug(bool status);

        static bool get_debug();

        static const RMIStats& get_stats();
    };
}

#endif // MADNESS_WORLD_WORLDRMI_H__INCLUDED
