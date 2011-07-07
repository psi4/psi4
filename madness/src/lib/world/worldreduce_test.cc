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

  $Id $
*/

#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <world/worldreduce.h>
#include <world/functional.h>
#include <world/functional.h>
#include <world/deferred_deleter.h>
#include <gtest/gtest.h>

madness::World* pworld;

using madness::WorldReduce;
using madness::detail::ReductionInterface;
using madness::detail::GroupReduction;

// WorldReduce must be used in the same way WorldObject is used.
// Note: This object must be stored in a shared pointer and the deleter must be
// a DeferredDeleter.
class Reducer : public madness::WorldReduce<Reducer, std::size_t> {
    typedef WorldReduce<Reducer, std::size_t> WorldReducer_;
public:

    Reducer(madness::World& w) :
        WorldReducer_(w)
    { process_pending(); }

    virtual ~Reducer() { }
};

template <typename T>
T add(const T& t1, const T& t2) { return t1 + t2; }

namespace {

    class WorldReduceTest : public ::testing::Test {
    public:
        WorldReduceTest() : reducer(new Reducer(*pworld), madness::DeferredDeleter<Reducer>(*pworld)) {
            for(ProcessID r = 0; r < pworld->size(); ++r) {
                if((r % 2) == 0)
                    even.push_back(r);

                all.push_back(r);
            }
        }

        virtual ~WorldReduceTest() { }

        ProcessID parent(const ProcessID node, const std::vector<ProcessID>& grp, const ProcessID root, unsigned int& count) {
            if(node == -1)
                return -1;

            ++count;
            GroupReduction<int> gr;
            gr.set_group(node, grp.begin(), grp.end(), root);

            unsigned int gr_count = 1;

            // Check that the children know who their parents are.
            if(gr.child0() != -1) {
                EXPECT_EQ(node, parent(gr.child0(), grp, root, count));
                ++gr_count;
            }
            if(gr.child1() != -1) {
                EXPECT_EQ(node, parent(gr.child1(), grp, root, count));
                ++gr_count;
            }

            // Check that the node count is correct.
            // 1 = no children
            // 2 = 1 child
            // 3 = 2 children
            EXPECT_EQ(gr_count, gr.count());

            if(gr.parent() == -1) {
                // test that the tree size is correct.
                EXPECT_EQ(gr.size(), count);

                // Check that root is correctly reported.
                EXPECT_TRUE(gr.is_root());
            } else {
                // Check that root is correctly reported.
                EXPECT_FALSE(gr.is_root());
            }

            return gr.parent();
        }

        int count_node(const ProcessID node, const std::vector<ProcessID>& grp, const ProcessID root, const ProcessID check_node) {
            if(node == -1)
                return 0;

            GroupReduction<int> gr;
            gr.set_group(node, grp.begin(), grp.end(), root);

            return (node == check_node ? 1 : 0)
                    + count_node(gr.child0(), grp, root, check_node)
                    + count_node(gr.child1(), grp, root, check_node);
        }

        std::vector<ProcessID> all;
        std::vector<ProcessID> even;
        std::shared_ptr<Reducer> reducer;
    };

    TEST_F(WorldReduceTest, Construct) {
        EXPECT_NO_THROW( GroupReduction<int> r );
        EXPECT_NO_THROW( Reducer x(*pworld) );
    }


    TEST_F(WorldReduceTest, Binary_Tree) {

        std::vector<ProcessID> grp;

        // Check group sizes from 0 to 20.
        for(int i = 0; i < 20; ++i) {
            grp.push_back(i);

            // Check for root at each node
            for(ProcessID root = 0; root < i; ++root) {
                unsigned int count = 0;

                // Check that children know who their parents are.
                EXPECT_EQ(-1, parent(root, grp, root, count));

                // Check that count equals group size.
                EXPECT_EQ(grp.size(),count);

                // Check that each node only appears once.
                for(std::vector<ProcessID>::const_iterator it = grp.begin(); it != grp.end(); ++it)
                    EXPECT_EQ(1,count_node(root, grp, root, *it));
            }

        }

        grp.resize(0);

        // Check group sizes from 0 to 40 for even nodes.
        for(int i = 0; i < 40; i += 2) {
            grp.push_back(i);

            // Check for root at each node
            for(ProcessID root = 0; root < i; root += 2) {
                unsigned int count = 0;

                // Check that children know who their parents are.
                EXPECT_EQ(-1, parent(root, grp, root, count));

                // Check that the number of nodes in the tree is equal to the group size.
                EXPECT_EQ(grp.size(),count);

                // Check that each node only appears once.
                for(std::vector<ProcessID>::const_iterator it = grp.begin(); it != grp.end(); ++it)
                    EXPECT_EQ(1, count_node(root, grp, root, *it));
            }

        }

    }

    TEST_F(WorldReduceTest, Reduce_All) {
        // Setup the reduction group
        // reduce(
        //      group key,
        //      local reduction value,
        //      reduction operation,
        //      first element in the reduction group list,
        //      last element in the reduction group list,
        //      root node)
        madness::Future<ProcessID> result = reducer->reduce(0, pworld->rank(),
                & add<ProcessID>, all.begin(), all.end(), 0);

        // The final value is available only on the root node.
        if(pworld->rank() == 0) {
            ProcessID sum = 0;
            for(std::vector<ProcessID>::const_iterator it = all.begin(); it != all.end(); ++it)
                sum += *it;
            EXPECT_EQ(sum, result.get());
        }

    }

    TEST_F(WorldReduceTest, Reduce_Even) {
        // You can only setup the reduction on nodes that are included in the group.
        if((pworld->rank() % 2) == 0) {
            madness::Future<ProcessID> result = reducer->reduce(1, pworld->rank(),
                    & add<ProcessID>, even.begin(), even.end(), 0);

            if(pworld->rank() == 0) {
                ProcessID sum = 0;
                for(std::vector<ProcessID>::const_iterator it = even.begin(); it != even.end(); ++it)
                    sum += *it;
                EXPECT_EQ(sum, result.get());
            }
        }

    }

    TEST_F(WorldReduceTest, Reduce_All_With_Future) {
        madness::Future<ProcessID> local_data;
        madness::Future<ProcessID> result = reducer->reduce(0, local_data,
                & add<ProcessID>, all.begin(), all.end(), 0);

        // Check that the result has not been set
        EXPECT_FALSE(result.probe());

        // set the data and fence
        local_data.set(pworld->rank());
        pworld->gop.fence();

        // Check that the reduction is complete on this node
        EXPECT_TRUE(result.probe());

        // The final value is available only on the root node.
        if(pworld->rank() == 0) {
            ProcessID sum = 0;
            for(std::vector<ProcessID>::const_iterator it = all.begin(); it != all.end(); ++it)
                sum += *it;
            EXPECT_EQ(sum, result.get());
        }

    }

    TEST_F(WorldReduceTest, Reduce_Even_With_Futures) {
        // You can only setup the reduction on nodes that are included in the group.
        if((pworld->rank() % 2) == 0) {
            madness::Future<ProcessID> local_data;
            madness::Future<ProcessID> result = reducer->reduce(1, local_data,
                    & add<ProcessID>, even.begin(), even.end(), 0);

            // Check that the result has not been set
            EXPECT_FALSE(result.probe());

            // set the data and wait for results
            local_data.set(pworld->rank());
            result.get();

            // Check that the reduction is complete on this node
            EXPECT_TRUE(result.probe());

            // The final value is available only on the root node.

            if(pworld->rank() == 0) {
                ProcessID sum = 0;
                for(std::vector<ProcessID>::const_iterator it = even.begin(); it != even.end(); ++it)
                    sum += *it;
                EXPECT_EQ(sum, result.get());
            }
        }

    }
}

int main(int argc, char **argv) {
    madness::initialize(argc,argv);
    int status = 0;
    {
        madness::World world(MPI::COMM_WORLD);
        pworld = &world;

        if (world.rank()) madness::redirectio(world);
        world.args(argc,argv);
        world.gop.fence();

        ::testing::InitGoogleTest(&argc, argv);
        status = RUN_ALL_TESTS();

        world.gop.fence();
    }
    madness::finalize();
    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the world test code\n";
    return 0;
}

#endif
