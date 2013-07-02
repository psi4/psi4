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

  $Id: safempi.h 2370 2011-06-15 00:21:35Z rjharrison $
*/
#ifndef MADNESS_WORLD_SAFEMPI_H__INCLUDED
#define MADNESS_WORLD_SAFEMPI_H__INCLUDED

/// \file safempi.h
/// \brief Serializes calls to MPI in case it does not support THREAD_MULTIPLE

// #ifdef STUBOUTMPI
// #include <world/stubmpi.h>
// #else
#include <mpi.h>
// #endif

#include <world/worldmutex.h>
#include <world/typestuff.h>
#include <world/enable_if.h>

namespace SafeMPI {
    
#ifndef HAVE_IBMBGP
#define SERIALIZE_MPI
#endif
    
#ifdef SERIALIZE_MPI
    extern madness::SCALABLE_MUTEX_TYPE charon;      // Inside safempi.cc
#define SAFE_MPI_GLOBAL_MUTEX madness::ScopedMutex<madness::SCALABLE_MUTEX_TYPE> obolus(SafeMPI::charon)
#else
#define SAFE_MPI_GLOBAL_MUTEX
#endif
    
    /// tags in [1,999] ... allocated once by unique_reserved_tag
    ///
    /// tags in [1000,1023] ... statically assigned here
    ///
    /// tags in [1024,4095] ... allocated round-robin by unique_tag
    ///
    /// tags in [4096,MPI::TAG_UB] ... not used/managed by madness
    
    static const int RMI_TAG = 1023;
    static const int RMI_HUGE_ACK_TAG = 1022;
    static const int RMI_HUGE_DAT_TAG = 1021;
    static const int MPIAR_TAG = 1001;
    static const int DEFAULT_SEND_RECV_TAG = 1000;
    
    class Request : private MPI::Request {
    public:
        Request() : MPI::Request() {}
        
        Request(const MPI::Request& request) : MPI::Request(request) {}
        
        bool Test() {
            SAFE_MPI_GLOBAL_MUTEX;
            return MPI::Request::Test();
        }
        
        bool Test_got_lock_already() {
            return MPI::Request::Test();
        }
        
        static bool Testany(int n, SafeMPI::Request* request, int& ind) {
            SAFE_MPI_GLOBAL_MUTEX;
            return MPI::Request::Testany(n, static_cast<MPI::Request*>(request), ind);
        }
        
        static int Testsome(int n, SafeMPI::Request* request, int* ind, MPI::Status* status) {
            SAFE_MPI_GLOBAL_MUTEX;
            return MPI::Request::Testsome(n, static_cast<MPI::Request*>(request), ind, status);
        }
        
        static int Testsome(int n, SafeMPI::Request* request, int* ind) {
            SAFE_MPI_GLOBAL_MUTEX;
            return MPI::Request::Testsome(n, static_cast<MPI::Request*>(request), ind);
        }
        
        // static int Testsome(int n, MPI::Request* request, int* ind, MPI::Status* status) {
        //     SAFE_MPI_GLOBAL_MUTEX;
        //     return MPI::Request::Testsome(n, request, ind, status);
        // }
    };
    
    class Intracomm {
        MPI::Intracomm& comm;
        int me;
        int numproc;
        
    public:
        Intracomm(MPI::Intracomm& c) : comm(c) {
            SAFE_MPI_GLOBAL_MUTEX;
            me = comm.Get_rank();
            numproc = comm.Get_size();
        }
        
        MPI::Intracomm Create(const MPI::Group& group) const {
            SAFE_MPI_GLOBAL_MUTEX;
            return comm.Create(group);
        }
        
        MPI::Group Get_group() const {
            SAFE_MPI_GLOBAL_MUTEX;
            return comm.Get_group();
        }

        int Get_rank() const { return me; }

        int Get_size() const { return numproc; }
        
        ::SafeMPI::Request Isend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            SAFE_MPI_GLOBAL_MUTEX;
            return comm.Isend(buf,count,datatype,dest,tag);
        }
        
        ::SafeMPI::Request Irecv(void* buf, size_t count, const MPI::Datatype& datatype, int src, int tag) const {
            SAFE_MPI_GLOBAL_MUTEX;
            return comm.Irecv(buf, count, datatype, src, tag);
        }
        
        void Send(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Send(buf,count,datatype,dest,tag);
        }
        
        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag, MPI::Status& status) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Recv(buf,count,datatype,source,tag,status);
        }
        
        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Recv(buf,count,datatype,source,tag);
        }
        
        void Bcast(void* buf, size_t count, const MPI::Datatype& datatype, int root) const {
            SAFE_MPI_GLOBAL_MUTEX;
            return comm.Bcast(buf, count, datatype, root);
        }
        
        void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op, int root) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Reduce(sendbuf, recvbuf, count, datatype, op, root);
        }
        
        void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Allreduce(sendbuf, recvbuf, count, datatype, op);
        }
        void Get_attr(int key, void* value) const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Get_attr(key, value);
        }
        
        void Abort(int code=1) const {
            comm.Abort(code);
        }
        
        void Barrier() const {
            SAFE_MPI_GLOBAL_MUTEX;
            comm.Barrier();
        }
        
        /// Returns a unique tag for temporary use (1023<tag<=4095)
        
        /// These tags are intended for one time use to avoid tag
        /// collisions with other messages around the same time period.
        /// It simply increments/wraps a counter and returns the next
        /// legal value.
        ///
        /// So that send and receiver agree on the tag all processes
        /// need to call this routine in the same sequence.
        static int unique_tag() {
            SAFE_MPI_GLOBAL_MUTEX;
            static volatile int tag = 1024;
            int result = tag++;
            if (tag >= 4095) tag = 1024;
            return result;
        }
        
        /// Returns a unique tag reserved for long-term use (0<tag<1000)
        
        /// Get a tag from this routine for long-term/repeated use.
        ///
        /// Tags in [1000,1023] are statically assigned.
        static int unique_reserved_tag() {
            SAFE_MPI_GLOBAL_MUTEX;
            static volatile int tag = 1;
            int result = tag++;
            if (result >= 1000) MADNESS_EXCEPTION( "too many reserved tags in use" , result );
            return result;
        }
        
        // Below here are extensions to MPI but have to be here since overload resolution
        // is not applied across different class scopes ... hence even having Isend with
        // a different signature in a derived class (WorldMpiInterface) hides this interface.
        // Thus, they must all be here.
        //
        // !! All of the routines below call the protected interfaces provided above.
        // !! Please ensure any additional routines follow this convention.
        /// Isend one element ... disabled for pointers to reduce accidental misuse.
        template <class T>
        typename madness::enable_if_c< !std::is_pointer<T>::value, SafeMPI::Request>::type
        Isend(const T& datum, int dest, int tag=DEFAULT_SEND_RECV_TAG) const {
            return Isend(&datum, sizeof(T), MPI::BYTE, dest, tag);
        }
        
        /// Async receive data of up to lenbuf elements from process dest
        template <class T>
        SafeMPI::Request
        Irecv(T* buf, int count, int source, int tag=DEFAULT_SEND_RECV_TAG) const {
            return Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
        }
        
        
        /// Async receive datum from process dest with default tag=1
        template <class T>
        typename madness::enable_if_c< !std::is_pointer<T>::value, SafeMPI::Request>::type
        Irecv(T& buf, int source, int tag=DEFAULT_SEND_RECV_TAG) const {
            return Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
        }
        
        
        /// Send array of lenbuf elements to process dest
        template <class T>
        void Send(const T* buf, long lenbuf, int dest, int tag=DEFAULT_SEND_RECV_TAG) const {
            Send((void*)buf, lenbuf*sizeof(T), MPI::BYTE, dest, tag);
        }
        
        
        /// Send element to process dest with default tag=1001
        
        /// Disabled for pointers to reduce accidental misuse.
        template <class T>
        typename madness::enable_if_c< !std::is_pointer<T>::value, void>::type
        Send(const T& datum, int dest, int tag=DEFAULT_SEND_RECV_TAG) const {
            Send((void*)&datum, sizeof(T), MPI::BYTE, dest, tag);
        }
        
        
        /// Receive data of up to lenbuf elements from process dest
        template <class T>
        void
        Recv(T* buf, long lenbuf, int src, int tag) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag);
        }
        
        /// Receive data of up to lenbuf elements from process dest with status
        template <class T>
        void
        Recv(T* buf, long lenbuf, int src, int tag, MPI::Status& status) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag, status);
        }
        
        
        /// Receive datum from process src
        template <class T>
        typename madness::enable_if_c< !std::is_pointer<T>::value, void>::type
        Recv(T& buf, int src, int tag=DEFAULT_SEND_RECV_TAG) const {
            Recv(&buf, sizeof(T), MPI::BYTE, src, tag);
        }
        
        
        /// MPI broadcast an array of count elements
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
        void Bcast(T* buffer, int count, int root) const {
            Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
        }
        
        
        /// MPI broadcast a datum
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
        typename madness::enable_if_c< !std::is_pointer<T>::value, void>::type
        Bcast(T& buffer, int root) const {
            Bcast(&buffer, sizeof(T), MPI::BYTE,root);
        }
        
        int rank() const { return Get_rank(); }
        
        int nproc() const { return Get_size(); }
        
        int size() const { return Get_size(); }
        
        
        /// Construct info about a binary tree with given root
        
        /// Constructs a binary tree spanning the communicator with
        /// process root as the root of the tree.  Returns the logical
        /// parent and children in the tree of the calling process.  If
        /// there is no parent/child the value -1 will be set.
        void binary_tree_info(int root, int& parent, int& child0, int& child1);
        
    };
}

#endif // MADNESS_WORLD_SAFEMPI_H__INCLUDED
