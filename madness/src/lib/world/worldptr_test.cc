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

#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#define MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE 0
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/worldptr.h>
#include <world/world.h>
#include <world/worldobj.h>
#include <world/bufar.h>

#include <gtest/gtest.h>

madness::World* pworld;

using madness::detail::WorldPtr;

namespace {
    class WorldPtrTest : public ::testing::Test {
    public:
        WorldPtrTest() : p(*pworld, i.get()), p0() {
            // You can do set-up work for each test here.
        }

        virtual ~WorldPtrTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case for World.
        static const std::shared_ptr<int> i;
        WorldPtr<int> p;
        WorldPtr<int> p0;

        class Object {
        public:
            int i;
        };

        template <typename T>
        class XferPtr : public madness::WorldObject<XferPtr<T> > {
            madness::Void set_ptr(const WorldPtr<T>& p, const madness::uniqueidT & id, bool away) {
                if(away)
                    remote_ptr.set(p);
                else
                    return_ptr.set(p);
                return madness::None;
            }
        public:
            madness::Future<WorldPtr<T> > remote_ptr;
            madness::Future<WorldPtr<T> > return_ptr;

            explicit XferPtr(madness::World& w) :
                    madness::WorldObject<XferPtr<T> >(w)
            {
                madness::WorldObject<XferPtr<T> >::process_pending();
            }

            void xfer(const ProcessID dest, const WorldPtr<T>& p, bool away) const {
                madness::WorldObject<XferPtr<T> >::send(dest, & XferPtr<T>::set_ptr,
                        p, madness::WorldObject<XferPtr<T> >::id(), away);
            }
        };
    };

    const std::shared_ptr<int> WorldPtrTest::i(new int(1));

    TEST_F(WorldPtrTest, DefaultConstructor) {
        // Check for correct initialization
        EXPECT_EQ(p0.owner(), -1);
        EXPECT_EQ(p0.get_worldid(), std::numeric_limits<unsigned long>::max());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p0.get(), madness::MadnessException);
        EXPECT_THROW(p0.get_world(), madness::MadnessException);
#endif
        EXPECT_FALSE(p0);
        EXPECT_TRUE(! p0);
        EXPECT_FALSE(p0.has_owner());
        EXPECT_FALSE(p0.is_local());

    }

    TEST_F(WorldPtrTest, PrimaryConstructor) {
        // Check for correct initialization
        EXPECT_EQ(p.get(), i.get());
        EXPECT_EQ(p.owner(), pworld->rank());
        EXPECT_EQ(p.get_worldid(), pworld->id());
        EXPECT_EQ(&(p.get_world()), pworld);
        EXPECT_TRUE(p);
        EXPECT_FALSE(! p);
        EXPECT_TRUE(p.has_owner());
        EXPECT_TRUE(p.is_local());
    }

    TEST_F(WorldPtrTest, CopyConstructor) {
        WorldPtr<int> c(p);

        // Check for correct initialization
        EXPECT_EQ(c.get(), p.get());
        EXPECT_EQ(c.owner(), p.owner());
        EXPECT_EQ(c.get_worldid(), pworld->id());
        EXPECT_EQ(&(c.get_world()), pworld);
        EXPECT_TRUE(c);
        EXPECT_FALSE(! c);
        EXPECT_TRUE(c.has_owner());
        EXPECT_TRUE(c.is_local());
    }

    TEST_F(WorldPtrTest, CopyConversionConstructor) {
        WorldPtr<const int> c(p);

        // Check for correct initialization
        EXPECT_EQ(c.get(), p.get());
        EXPECT_EQ(c.owner(), p.owner());
        EXPECT_EQ(c.get_worldid(), pworld->id());
        EXPECT_EQ(&(c.get_world()), pworld);
        EXPECT_TRUE(c);
        EXPECT_FALSE(! c);
        EXPECT_TRUE(c.has_owner());
        EXPECT_TRUE(c.is_local());
    }

    TEST_F(WorldPtrTest, Assignment) {
        WorldPtr<int> c;

        c = p;

        // Check for correct initialization
        EXPECT_EQ(c.get(), p.get());
        EXPECT_EQ(c.owner(), p.owner());
        EXPECT_EQ(c.get_world().id(), p.get_world().id());
    }

    TEST_F(WorldPtrTest, ConversionAssignment) {
        WorldPtr<const int> c;

        c = p;

        // Check for correct initialization
        EXPECT_EQ(c.get(), p.get());
        EXPECT_EQ(c.owner(), p.owner());
        EXPECT_EQ(c.get_world().id(), p.get_world().id());
    }

    TEST_F(WorldPtrTest, PointerAccess) {

        Object obj;
        obj.i = 1;
        WorldPtr<Object> p(*pworld, &obj);

        // check for correct pointer and member access
        EXPECT_EQ(p.get(), &obj);
        EXPECT_EQ((*p).i, obj.i);
        EXPECT_EQ(&(*p), &obj);
        EXPECT_EQ(p->i, obj.i);
    }

    TEST_F(WorldPtrTest, LocalOwnership) {
        // Check type conversion to bool
        EXPECT_TRUE(p.has_owner());
        EXPECT_FALSE(p0.has_owner());

        EXPECT_TRUE(p.is_local());
        EXPECT_FALSE(p0.is_local());

        EXPECT_EQ(p.owner(), pworld->rank());
        EXPECT_EQ(p0.owner(), -1);
    }

    TEST_F(WorldPtrTest, WorldAccess) {
        // Check access to world
        EXPECT_EQ(&(p.get_world()), pworld);

        // Check for assertion
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p0.get_world(), madness::MadnessException);
#endif

        // Check for correct world id
        EXPECT_EQ(p0.get_worldid(), std::numeric_limits<unsigned long>::max());
        EXPECT_EQ(p.get_worldid(), pworld->id());
    }

    TEST_F(WorldPtrTest, BooleanOperations) {

        WorldPtr<int> c(p);
        WorldPtr<int> c0(p0);

        // Check type conversion to bool
        bool bp = p;
        bool bp0 = p0;
        EXPECT_TRUE(bp);
        EXPECT_FALSE(bp0);

        bool bnp = ! p;
        bool bnp0 = ! p0;
        EXPECT_FALSE(bnp);
        EXPECT_TRUE(bnp0);

        // Check equality operator
        EXPECT_FALSE(p == p0);
        EXPECT_FALSE(p0 == p);

        EXPECT_TRUE(p == c);
        EXPECT_TRUE(p0 == c0);

        // Check inequality operator
        EXPECT_TRUE(p != p0);
        EXPECT_TRUE(p0 != p);

        EXPECT_FALSE(p != c);
        EXPECT_FALSE(p0 != c0);

        // Check less-than operator
        EXPECT_FALSE(p < p0);
        EXPECT_TRUE(p0 < p);

        EXPECT_FALSE(p < c);
        EXPECT_FALSE(p0 < c0);

        // Check greater-than operator
        EXPECT_TRUE(p > p0);
        EXPECT_FALSE(p0 > p);

        EXPECT_FALSE(p > c);
        EXPECT_FALSE(p0 > c0);

        // Check less-than-or-equal-to operator
        EXPECT_FALSE(p <= p0);
        EXPECT_TRUE(p0 <= p);

        EXPECT_TRUE(p <= c);
        EXPECT_TRUE(p0 <= c0);

        // Check greater-than-or-equal-to operator
        EXPECT_TRUE(p >= p0);
        EXPECT_FALSE(p0 >= p);

        EXPECT_TRUE(p >= c);
        EXPECT_TRUE(p0 >= c0);
    }


    TEST_F(WorldPtrTest, Swap) {
        // Verify initial conditions
        EXPECT_EQ(p.get(), i.get());
        EXPECT_EQ(p.owner(), pworld->rank());
        EXPECT_EQ(&(p.get_world()), pworld);
        EXPECT_EQ(p.get_worldid(), pworld->id());

        EXPECT_EQ(p0.owner(), -1);
        EXPECT_EQ(p0.get_worldid(), std::numeric_limits<unsigned long>::max());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p0.get(), madness::MadnessException);
        EXPECT_THROW(p0.get_world(), madness::MadnessException);
#endif

        p.swap(p0);

        // Check member swap post conditions
        EXPECT_EQ(p0.get(), i.get());
        EXPECT_EQ(p0.owner(), pworld->rank());
        EXPECT_EQ(&(p0.get_world()), pworld);
        EXPECT_EQ(p0.get_worldid(), pworld->id());

        EXPECT_EQ(p.owner(), -1);
        EXPECT_EQ(p.get_worldid(), std::numeric_limits<unsigned long>::max());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p.get(), madness::MadnessException);
        EXPECT_THROW(p.get_world(), madness::MadnessException);
#endif

        madness::detail::swap(p,p0);

        // Check free function swap post conditions conditions
        EXPECT_EQ(p.get(), i.get());
        EXPECT_EQ(p.owner(), pworld->rank());
        EXPECT_EQ(&(p.get_world()), pworld);
        EXPECT_EQ(p.get_worldid(), pworld->id());

        EXPECT_EQ(p0.owner(), -1);
        EXPECT_EQ(p0.get_worldid(), std::numeric_limits<unsigned long>::max());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p0.get(), madness::MadnessException);
        EXPECT_THROW(p0.get_world(), madness::MadnessException);
#endif

        std::swap(p,p0);

        // Check std::swap post conditions
        EXPECT_EQ(p0.get(), i.get());
        EXPECT_EQ(p0.owner(), pworld->rank());
        EXPECT_EQ(&(p0.get_world()), pworld);
        EXPECT_EQ(p0.get_worldid(), pworld->id());

        EXPECT_EQ(p.owner(), -1);
        EXPECT_EQ(p.get_worldid(), std::numeric_limits<unsigned long>::max());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(p.get(), madness::MadnessException);
        EXPECT_THROW(p.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldPtrTest, Serialize) {
        // Serialize 2 pointers to a buffer
        unsigned char buf[10*sizeof(WorldPtr<int>)];
        madness::archive::BufferOutputArchive oar(buf,sizeof(buf));
        oar & p & p0;
        std::size_t nbyte = oar.size();
        oar.close();

        // Deserialize 2 pointers from a buffer
        WorldPtr<int> a;
        WorldPtr<int> a0(p);
        madness::archive::BufferInputArchive iar(buf,nbyte);
        iar & a & a0;
        iar.close();

        // check that serialized and deserialized pointers match the original.
        EXPECT_EQ(a.get(), p.get());
        EXPECT_EQ(a.owner(), p.owner());
        EXPECT_EQ(&(a.get_world()), &(p.get_world()));
        EXPECT_EQ(a.get_worldid(), p.get_worldid());

        EXPECT_TRUE(a0 == p0);
        EXPECT_EQ(a0.owner(), p0.owner());
        EXPECT_EQ(a0.get_worldid(), p0.get_worldid());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(a0.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldPtrTest, RemotePtr) {

        if(pworld->size() > 1) {
            ProcessID send = (pworld->rank() != (pworld->size() - 1) ? pworld->rank() + 1 : 0);
            ProcessID recv = (pworld->rank() != 0 ? pworld->rank() - 1 : (pworld->size() - 1) );

            // Send the pointer to the next process
            XferPtr<int> xfer_wobj(*pworld);
            xfer_wobj.xfer(send, p, true);

            // wait for the remote pointer
            WorldPtr<int> remote = xfer_wobj.remote_ptr.get();

            EXPECT_FALSE(remote == p);

            // Check for locality
            EXPECT_FALSE(remote.is_local());

            // Check ownership
            EXPECT_TRUE(remote.has_owner());
            EXPECT_EQ(remote.owner(), recv);

            // Check world and world id of remote pointer
            EXPECT_EQ(&(remote.get_world()), pworld);
            EXPECT_EQ(remote.get_worldid(), pworld->id());

            // Check boolean conversions
            EXPECT_TRUE(remote);
            EXPECT_FALSE(! remote);

#ifdef MADNESS_ASSERTIONS_THROW
            EXPECT_THROW(remote.get(), madness::MadnessException);
            EXPECT_THROW(*remote, madness::MadnessException);
            // Should check arrow operator too but the assertions are the same as
            // operator* so we should get the same throw conditions
#endif


            xfer_wobj.xfer(recv, remote, false);

            // wait for the original to return.
            WorldPtr<int> back = xfer_wobj.return_ptr.get();

            EXPECT_TRUE(back == p);

            // Check for locality
            EXPECT_TRUE(back.is_local());

            // Check ownership
            EXPECT_TRUE(back.has_owner());
            EXPECT_EQ(back.owner(), pworld->rank());

            // Check world and world id of remote pointer
            EXPECT_EQ(&(back.get_world()), pworld);
            EXPECT_EQ(back.get_worldid(), pworld->id());

            // Check boolean conversions
            EXPECT_TRUE(back);
            EXPECT_FALSE(! back);


            EXPECT_EQ(back.get(), p.get());
            EXPECT_EQ(*back, *p);
        }
    }

    TEST_F(WorldPtrTest, RemoteNullPtr) {

        if(pworld->size() > 1) {

            ProcessID send = (pworld->rank() != (pworld->size() - 1) ? pworld->rank() + 1 : 0);
            ProcessID recv = (pworld->rank() != 0 ? pworld->rank() - 1 : (pworld->size() - 1) );

            // Send the pointer to the next process
            XferPtr<int> xfer_wobj(*pworld);
            xfer_wobj.xfer(send, p0, true);

            // wait for the remote pointer
            WorldPtr<int> remote = xfer_wobj.remote_ptr.get();

            EXPECT_TRUE(remote == p0);

            // Check for locality
            EXPECT_FALSE(remote.is_local());

            // Check ownership
            EXPECT_FALSE(remote.has_owner());
            EXPECT_EQ(remote.owner(), -1);

            // Check world and world id of remote pointer

            EXPECT_EQ(remote.get_worldid(), std::numeric_limits<unsigned long>::max());

            // Check boolean conversions
            EXPECT_FALSE(remote);
            EXPECT_TRUE(! remote);

#ifdef MADNESS_ASSERTIONS_THROW
            EXPECT_THROW(remote.get(), madness::MadnessException);
            EXPECT_THROW(remote.get_world(), madness::MadnessException);
            EXPECT_THROW(*remote, madness::MadnessException);
            // Should check arrow operator too but the assertions are the same as
            // operator* so we should get the same throw conditions
#endif


            xfer_wobj.xfer(recv, remote, false);

            // wait for the original to return.
            WorldPtr<int> back = xfer_wobj.return_ptr.get();

            EXPECT_TRUE(back == p0);

            // Check for locality
            EXPECT_FALSE(back.is_local());

            // Check ownership
            EXPECT_FALSE(back.has_owner());
            EXPECT_EQ(back.owner(), -1);

            // Check world and world id of remote pointer
#ifdef MADNESS_ASSERTIONS_THROW
            EXPECT_THROW(back.get(), madness::MadnessException);
            EXPECT_THROW(back.get_world(), madness::MadnessException);
            EXPECT_THROW(*back, madness::MadnessException);
#endif
            EXPECT_EQ(back.get_worldid(), std::numeric_limits<unsigned long>::max());

            // Check boolean conversions
            EXPECT_FALSE(back);
            EXPECT_TRUE(! back);
        }
    }


} // namespace

int main(int argc, char **argv) {
    madness::initialize(argc,argv);
    madness::World world(MPI::COMM_WORLD);
    pworld = &world;

    if (world.rank()) madness::redirectio(world);
    world.args(argc,argv);
    world.gop.fence();

    ::testing::InitGoogleTest(&argc, argv);
    int status = RUN_ALL_TESTS();

    world.gop.fence();

    madness::finalize();
    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the WorldPtr test code\n";
    return 0;
}

#endif
