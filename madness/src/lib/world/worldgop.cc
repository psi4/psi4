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


  $Id: $
*/

#include <world/worldgop.h>
#include <world/world.h> // for World, WorldTaskQueue, and WorldAmInterface

namespace madness {



    void WorldGopInterface::await(SafeMPI::Request& req) { World::await(req); }

    // In the World constructor can ONLY rely on MPI and MPI being initialized
    WorldGopInterface::WorldGopInterface(World& world)
            : mpi(world.mpi)
            , am(world.am)
            , taskq(world.taskq)
            , deferred(new detail::DeferredCleanup())
            , debug(false)
    { }

    WorldGopInterface::~WorldGopInterface() {
        deferred->destroy(true);
        deferred->do_cleanup();
    }

    /// Set debug flag to new value and return old value
    bool WorldGopInterface::set_debug(bool value) {
        bool status = debug;
        debug = value;
        return status;
    }

    /// Synchronizes all processes in communicator ... does NOT fence pending AM or tasks
    void WorldGopInterface::barrier() {
        long i = rank();
        sum(i);
        if (i != size()*(size()-1)/2) error("bad value after sum in barrier");
    }


    /// Synchronizes all processes in communicator AND globally ensures no pending AM or tasks

    /// Runs Dykstra-like termination algorithm on binary tree by
    /// locally ensuring ntask=0 and all am sent and processed,
    /// and then participating in a global sum of nsent and nrecv.
    /// Then globally checks that nsent=nrecv and that both are
    /// constant over two traversals.  We are then we are sure
    /// that all tasks and AM are processed and there no AM in
    /// flight.
    void WorldGopInterface::fence() {
        PROFILE_MEMBER_FUNC(WorldGopInterface);
        unsigned long nsent_prev=0, nrecv_prev=1; // invalid initial condition
        SafeMPI::Request req0, req1;
        ProcessID parent, child0, child1;
        mpi.binary_tree_info(0, parent, child0, child1);
        Tag gfence_tag = mpi.unique_tag();
        int npass = 0;

        //double start = wall_time();

        while (1) {
            uint64_t sum0[2]={0,0}, sum1[2]={0,0}, sum[2];
            if (child0 != -1) req0 = mpi.Irecv((void*) &sum0, sizeof(sum0), MPI::BYTE, child0, gfence_tag);
            if (child1 != -1) req1 = mpi.Irecv((void*) &sum1, sizeof(sum1), MPI::BYTE, child1, gfence_tag);
            taskq.fence();
            if (child0 != -1) World::await(req0);
            if (child1 != -1) World::await(req1);

            bool finished;
            uint64_t ntask1, nsent1, nrecv1, ntask2, nsent2, nrecv2;
            do {
                taskq.fence();

                // Since the number of outstanding tasks and number of AM sent/recv
                // don't share a critical section read each twice and ensure they
                // are unchanged to ensure that are consistent ... they don't have
                // to be current.

                ntask1 = taskq.size();
                nsent1 = am.nsent;
                nrecv1 = am.nrecv;

                __asm__ __volatile__ (" " : : : "memory");

                ntask2 = taskq.size();
                nsent2 = am.nsent;
                nrecv2 = am.nrecv;

                __asm__ __volatile__ (" " : : : "memory");

                finished = (ntask2==0) && (ntask1==0) && (nsent1==nsent2) && (nrecv1==nrecv2);
            }
            while (!finished);

            sum[0] = sum0[0] + sum1[0] + nsent2; // Must use values read above
            sum[1] = sum0[1] + sum1[1] + nrecv2;

            if (parent != -1) {
                req0 = mpi.Isend(&sum, sizeof(sum), MPI::BYTE, parent, gfence_tag);
                World::await(req0);
            }

            // While we are probably idle free unused communication buffers
            am.free_managed_buffers();

            //bool dowork = (npass==0) || (ThreadPool::size()==0);
            bool dowork = true;
            broadcast(&sum, sizeof(sum), 0, dowork);
            ++npass;

//            madness::print("GOPFENCE", npass, sum[0], nsent_prev, sum[1], nrecv_prev);

            if (sum[0]==sum[1] && sum[0]==nsent_prev && sum[1]==nrecv_prev) {
                break;
            }

//                 if (wall_time() - start > 1200.0) {
//                     std::cout << rank() << " FENCE " << nsent2 << " "
//                         << nsent_prev << " " << nrecv2 << " " << nrecv_prev
//                         << " " << sum[0] << " " << sum[1] << " " << npass
//                         << " " << taskq.size() << std::endl;
//                     std::cout.flush();
//                     //myusleep(1000);
//                     MADNESS_ASSERT(0);
//                 }

            nsent_prev = sum[0];
            nrecv_prev = sum[1];

        };
        am.free_managed_buffers(); // free up communication buffers
        deferred->do_cleanup();
    }


    /// Broadcasts bytes from process root while still processing AM & tasks

    /// Optimizations can be added for long messages
    void WorldGopInterface::broadcast(void* buf, size_t nbyte, ProcessID root, bool dowork) {
        SafeMPI::Request req0, req1;
        ProcessID parent, child0, child1;
        mpi.binary_tree_info(root, parent, child0, child1);
        Tag bcast_tag = mpi.unique_tag();

        //print("BCAST TAG", bcast_tag);

        if (parent != -1) {
            req0 = mpi.Irecv(buf, nbyte, MPI::BYTE, parent, bcast_tag);
            World::await(req0, dowork);
        }

        if (child0 != -1) req0 = mpi.Isend(buf, nbyte, MPI::BYTE, child0, bcast_tag);
        if (child1 != -1) req1 = mpi.Isend(buf, nbyte, MPI::BYTE, child1, bcast_tag);

        if (child0 != -1) World::await(req0, dowork);
        if (child1 != -1) World::await(req1, dowork);
    }

} // namespace madness
