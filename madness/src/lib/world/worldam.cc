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

#include <world/worldam.h>
#include <world/world.h>
#include <world/worldmpi.h>

namespace madness {



    WorldAmInterface::WorldAmInterface(World& world)
            : world(world)
            , rank(world.mpi.Get_rank())
            , nproc(world.mpi.Get_size())
            , cur_msg(0)
            , nsent(0)
            , nrecv(0)
            , map_to_comm_world(nproc)
    {
        lock();
        for (int i=0; i<NSEND; ++i) managed_send_buf[i] = 0;

        std::vector<int> fred(nproc);
        for (int i=0; i<nproc; ++i) fred[i] = i;
        MPI::Group::Translate_ranks(world.mpi.comm().Get_group(), nproc, &fred[0],
                                    MPI::COMM_WORLD.Get_group(), &map_to_comm_world[0]);

        // for (int i=0; i<nproc; ++i) {
        //     std::cout << "map " << i << " " << map_to_comm_world[i] << std::endl;
        // }

        unlock();
    }

    WorldAmInterface::~WorldAmInterface() {
        for (int i=0; i<NSEND; ++i) {
            while (!send_req[i].Test()) {
                myusleep(100);
            }
            free_managed_send_buf(i);
        }
    }

    void AmArg::set_world(World* world) const {
        worldid = world->id();
    }

    World* AmArg::get_world() const {
        return World::world_from_id(worldid);
    }

    void WorldAmInterface::increment_worldam_nrecv(World* world) {
        world->am.nrecv++;
    }

} // namespace madness
