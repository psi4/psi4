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


  $Id: worldam.h 2321 2011-05-23 00:49:51Z rjharrison $
*/

#ifndef MADNESS_WORLD_WORLDAM_H__INCLUDED
#define MADNESS_WORLD_WORLDAM_H__INCLUDED

/// \file worldam.h
/// \brief Implements active message layer for World on top of RMI layer

#include <world/bufar.h>
#include <world/worldrmi.h>
#include <vector>
#include <cstddef>

namespace madness {

    /*
      The RMI layer just does transport and does not know about World
      or even necessarily about MPI.  It also has no buffering or
      virtualization of resources.  In particular, we must be careful
      about having too many outstanding sends and active message
      handlers must be careful about what they do in order to avoid
      deadlock --- especially problematic is a handler trying to send
      a message.

      The WorldAM class provides a World-aware RMI capability that
      limits the number of outstanding sends and can optionally manage
      buffers.  The issue of what handlers can do safely is handled by
      the integration of WorldAM and the task queues ... if you have
      an operation that might

      - send messages

      - take a long time

      - consume a lot of stack/heap (e.g., recursive algorithm)

      then the right thing to do is to send a task rather than
      an active message.
    */

    class World;
    template <class Derived> class WorldObject;

    class AmArg;
    /// Type of AM handler functions
    typedef void (*am_handlerT)(const AmArg&);

    /// World active message that extends an RMI message
    class AmArg {
    private:
        friend class WorldAmInterface;
        template <class Derived> friend class WorldObject;

        friend AmArg* alloc_am_arg(std::size_t nbyte);

        unsigned char header[RMI::HEADER_LEN]; // !!!!!!!!!  MUST BE FIRST !!!!!!!!!!
        std::size_t nbyte;              // Size of user payload
        mutable unsigned long worldid;  // Id of associated world
        mutable am_handlerT func;       // User function to call
        mutable ProcessID src;          // Rank of process sending the message
        mutable unsigned int flags;     // Misc. bit flags

        // On 32 bit machine AmArg is HEADER_LEN+4+4+4+4+4=84 bytes
        // On 64 bit machine AmArg is HEADER_LEN+8+8+8+4+4=96 bytes

        // No copy constructor or assignment
        AmArg(const AmArg&);
        AmArg& operator=(const AmArg&);

        void set_src(ProcessID source) const { src = source; }

        // This is not inline in order to keep World opaque.
        void set_world(World* world) const;

        void set_func(am_handlerT handler) const { func = handler; }

        void set_size(std::size_t numbyte) { nbyte = numbyte; }

        void set_pending() const { flags |= 0x1ul; }

        bool is_pending() const { return flags & 0x1ul; }
        
        void clear_flags() const { flags = 0; }
        
        am_handlerT get_func() const { return func; }

        archive::BufferInputArchive make_input_arch() const {
            return archive::BufferInputArchive(buf(),size());
        }

        archive::BufferOutputArchive make_output_arch() const {
            return archive::BufferOutputArchive(buf(),size());
        }
        
    public:
        AmArg() {}

        /// Returns a pointer to the user's payload (aligned in same way as AmArg)
        unsigned char* buf() const { return (unsigned char*)(this) + sizeof(AmArg); }

        /// Returns the size of the user's payload
        std::size_t size() const { return nbyte; }

        /// Used to deserialize arguments from incoming message
        template <typename T>
        archive::BufferInputArchive operator&(T& t) const {
            return make_input_arch() & t;
        }

        /// Used to serialize arguments into outgoing message
        template <typename T>
        archive::BufferOutputArchive operator&(const T& t) const {
            return make_output_arch() & t;
        }

        /// For incoming AM gives the source process
        ProcessID get_src() const { return src; }

        // This is not inline in order to keep World opaque.
        /// For incoming AM gives the associated world
        World* get_world() const;
    };


    /// Allocates a new AmArg with nbytes of user data ... delete with free_am_arg
    inline AmArg* alloc_am_arg(std::size_t nbyte) {
        std::size_t narg = 1 + (nbyte+sizeof(AmArg)-1)/sizeof(AmArg);
        AmArg *arg = new AmArg[narg];
        arg->set_size(nbyte);
        return arg;
    }


    inline AmArg* copy_am_arg(const AmArg& arg) {
        AmArg* r = alloc_am_arg(arg.size());
        memcpy(r, &arg, arg.size()+sizeof(AmArg));
        return r;
    }

    /// Frees an AmArg allocated with alloc_am_arg
    inline void free_am_arg(AmArg* arg) {
        delete [] arg;
    }


    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h,
                             const I& i, const J& j) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d,e,f,g,h,i,j));
        *arg & a & b & c & d & e & f & g & h & i & j;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h,
                             const I& i) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d,e,f,g,h,i));
        *arg & a & b & c & d & e & f & g & h & i;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d,e,f,g,h));
        *arg & a & b & c & d & e & f & g & h;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d,e,f,g));
        *arg & a & b & c & d & e & f & g;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E, typename F>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,e,d,f));
        *arg & a & b & c & d & e & f;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D, typename E>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d, const E& e) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d,e));
        *arg & a & b & c & d & e;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C, typename D>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c, const D& d) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c,d));
        *arg & a & b & c & d;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B, typename C>
    inline AmArg* new_am_arg(const A& a, const B& b, const C& c) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b,c));
        *arg & a & b & c;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A, typename B>
    inline AmArg* new_am_arg(const A& a, const B& b) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a,b));
        *arg & a & b;
        return arg;
    }

    /// Convenience template for serializing arguments into a new AmArg
    template <typename A>
    inline AmArg* new_am_arg(const A& a) {
        AmArg* arg = alloc_am_arg(archive::bufar_size(a));
        *arg & a;
        return arg;
    }

    /// Implements AM interface
    class WorldAmInterface : private SCALABLE_MUTEX_TYPE {
        friend class WorldGopInterface;
    public:
        static const int MSG_LEN = RMI::MAX_MSG_LEN - sizeof(AmArg); ///< Max length of user payload in message
    private:
#ifdef HAVE_CRAYXT
        static const int NSEND = 512; ///< Max no. of pending sends
#else
        static const int NSEND = 32;///< Max no. of pending sends
#endif

        // Multiple threads are making their way thru here ... must be careful
        // to ensure updates are atomic and consistent

        AmArg* volatile managed_send_buf[NSEND]; ///< Managed buffers
        RMI::Request send_req[NSEND];   ///< Tracks in-flight messages

        World& world;            ///< The world which contains this instance of WorldAmInterface
        const ProcessID rank;
        const int nproc;
        volatile int cur_msg; ///< Index of next buffer to attempt to use
        volatile unsigned long nsent;     ///< Counts no. of AM sent for purpose of termination detection
        volatile unsigned long nrecv;     ///< Counts no. of AM received for purpose of termination detection

        std::vector<int> map_to_comm_world; ///< Maps rank in current MPI communicator to MPI::COMM_WORLD

        void free_managed_send_buf(int i) {
            // WE ASSUME WE ARE INSIDE A CRITICAL SECTION WHEN IN HERE
            if (managed_send_buf[i]) {
                free_am_arg(managed_send_buf[i]);
                managed_send_buf[i] = 0;
            }
        }

        /// Private: Finds/waits for a free send request
    int get_free_send_request() {
        // WE ASSUME WE ARE INSIDE A CRITICAL SECTION WHEN IN HERE
//             // Sequentially loop looking for next free request.
//             while (!send_req[cur_msg].Test()) {
//                 cur_msg++;
//                 if (cur_msg >= NSEND) cur_msg = 0;
//                 myusleep(5);
//             }

        // Wait for oldest request to complete
        while (!send_req[cur_msg].Test()) {
            // If the oldest message has still not completed then there is likely
            // severe network or end-point congestion, so pause for 100us in a rather
            // abitrary attempt to decreate the injection rate.  The server thread
            // is still polling every 1us (which is required to suck data off the net
            // and by some engines to ensure progress on sends).
            myusleep(100);
        }

        free_managed_send_buf(cur_msg);
        int result = cur_msg;
        cur_msg++;
        if (cur_msg >= NSEND) cur_msg = 0;

        return result;
    }

        // Not inline in order to keep World opaque
        static void increment_worldam_nrecv(World* world);

        /// This handles all incoming RMI messages for all instances
        static void handler(void *buf, std::size_t nbyte) {
            // It will be singled threaded since only the RMI receiver
            // thread will invoke it ... however note that nrecv will
            // be read by the main thread during fence operations.
            AmArg* arg = (AmArg*)(buf);
            am_handlerT func = arg->get_func();
            World* world = arg->get_world();
            MADNESS_ASSERT(arg->size() + sizeof(AmArg) == nbyte);
            MADNESS_ASSERT(world);
            MADNESS_ASSERT(func);
            func(*arg);
            //world->am.nrecv++;  // Must be AFTER execution of the function
            increment_worldam_nrecv(world);  // Must be AFTER execution of the function
        }
        
        /// Sends a non-blocking active message
        RMI::Request isend(ProcessID dest, am_handlerT op, const AmArg* arg, int attr, bool managed) {
            arg->set_world(&world);
            arg->set_src(rank);
            arg->set_func(op);
            arg->clear_flags(); // Is this the right place for this?
            
            MADNESS_ASSERT(arg->get_world());
            MADNESS_ASSERT(arg->get_func());
            
            // Map dest from world's communicator to comm_world
            dest = map_to_comm_world[dest];
            
            lock();    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            nsent++;
            int i = get_free_send_request();
            send_req[i] = RMI::isend(arg, arg->size()+sizeof(AmArg), dest, handler, attr);
            if (managed) managed_send_buf[i] = (AmArg*)(arg);
            unlock();  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            return send_req[i];
        }
        

    public:
        WorldAmInterface(World& world);

        virtual ~WorldAmInterface();

        /// Currently a noop
        void fence() {}

        /// Sends an unmanaged non-blocking active message
//        RMI::Request isend(ProcessID dest, am_handlerT op, const AmArg* arg, int attr=RMI::ATTR_ORDERED);

        /// Sends a managed non-blocking active message
    void send(ProcessID dest, am_handlerT op, const AmArg* arg, int attr=RMI::ATTR_ORDERED) {
        isend(dest, op, arg, attr, true);
    }
        /// Frees as many send buffers as possible
    void free_managed_buffers() {
        int ind[NSEND];
        lock(); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        int n = SafeMPI::Request::Testsome(NSEND, send_req, ind);
        if (n != MPI_UNDEFINED) {
            for (int i=0; i<n; ++i) {
                free_managed_send_buf(ind[i]);
            }
        }
        unlock(); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

//    RMI::Request isend(ProcessID dest, am_handlerT op, const AmArg* arg, int attr) {
//        std::cerr << "ISEND_ING AM\n";
//        return isend(dest, op, arg, attr, false);
//    }
        
    };
}

#endif // MADNESS_WORLD_WORLDAM_H__INCLUDED
