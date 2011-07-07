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

#include <world/worldtask.h>
//#include <world/worldmpi.h>

namespace madness {

    bool TaskInterface::debug = false;

    void TaskInterface::run(const TaskThreadEnv& env) { // This is what thread pool will invoke
        MADNESS_ASSERT(world);
        MADNESS_ASSERT(completion);
        World* w = const_cast<World*>(world);
        if (debug) std::cerr << w->rank() << ": Task " << (void*) this << " is now running" << std::endl;
        run(*w, env);
        if (debug) std::cerr << w->rank() << ": Task " << (void*) this << " has completed" << std::endl;
    }

    WorldTaskQueue::WorldTaskQueue(World& world)
            : world(world)
            , me(world.rank()) {
        nregistered = 0;
    }

    bool WorldTaskQueue::Stealer::operator()(PoolTaskInterface** pt) {
        madness::print("IN STEAL");
        PoolTaskInterface* t = *pt;
        if (t->is_stealable()) {
            TaskInterface* task = dynamic_cast<TaskInterface*>(t);
            if (task) {
                if (task->get_world()->id() == q.world.id()) {
                    madness::print("Stealing", (void *) task);
                    v.push_back(task);
                    *pt = 0; // Indicates task has been stolen
                }
            }
        }
        return true;
    }

    std::vector<TaskInterface*> WorldTaskQueue::steal(int nsteal) {
        std::vector<TaskInterface*> v;
        Stealer xxx(*this, v, nsteal);
        ThreadPool::instance()->scan(xxx);
        return v;
    }

    bool WorldTaskQueue::ProbeAllDone::operator()() const {
        if (cpu_time()-start > 1200) {
            for (int loop = 0; loop<3; ++loop) {
                std::cout << "HUNG Q? " << tq->size() << " " << ThreadPool::queue_size() << std::endl;
                std::cout.flush();
                myusleep(1000000);
            }
            MADNESS_ASSERT(cpu_time()-start < 1200);
        }
        return (tq->size() == 0);
    }

    void WorldTaskQueue::fence() {
        ProbeAllDone tester(this);
        do {
            world.await(tester);
        }
        while (nregistered);
    }

}  // namespace madness
