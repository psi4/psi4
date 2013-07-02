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
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <world/worlddc.h>

using namespace madness;
using namespace std;

struct Key {
    int k;

    Key() : k(-1) {}

    Key(int k) : k(k) {}

    hashT hash() const {
        return k;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & k;
    }

    bool operator==(const Key& b) const {
        return k==b.k;
    }
};

ostream& operator<<(ostream&s, const Key& key) {
    s << "Key(" << key.k << ")";
    return s;
}

struct Node {
    int k;

    Node() : k(-1) {}

    Node(int k) : k(k) {}

    int get() const {
        return k;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & k;
    }
};

ostream& operator<<(ostream&s, const Node& node) {
    s << "Node(" << node.k << ")";
    return s;
}


void test0(World& world) {
    WorldContainer<Key,Node> c(world);

    Key key1(1);
    Node node1(1);

    if (c.owner(key1) == world.rank()) c.replace(key1,node1);

    world.gop.fence();

    for (int i=0; i<10000; ++i)
        MADNESS_ASSERT(c.find(key1).get()->second.get() == 1);

    for (int i=3; i<100; ++i)
        MADNESS_ASSERT(c.find(Key(i)).get() == c.end());


    world.gop.fence();
}

class TestPmap : public WorldDCPmapInterface<int> {
private:
    const int nproc;
    const int shift;
public:
    TestPmap(World& world, int shift)
        : nproc(world.mpi.nproc())
        , shift(shift)
    { }

    ProcessID owner(const int& key) const {
        if (nproc == 1) return 0;
        return (key + shift)%nproc;
    }
};

void test1(World& world) {
    std::shared_ptr< WorldDCPmapInterface<int> > pmap0(new TestPmap(world, 0));
    std::shared_ptr< WorldDCPmapInterface<int> > pmap1(new TestPmap(world, 1));

    WorldContainer<int,double> c(world,pmap0), d(world,pmap0), e(world,pmap0);

    if (world.rank() == 0) {
        for (int i=0; i<100; ++i) {
            c.replace(i,i+1.0);
            d.replace(i,i+2.0);
            e.replace(i,i+3.0);
        }
    }

    pmap0->redistribute(world, pmap1);

    for (int i=0; i<100; ++i) {
        MADNESS_ASSERT(c.find(i).get()->second == (i+1.0));
        MADNESS_ASSERT(d.find(i).get()->second == (i+2.0));
        MADNESS_ASSERT(e.find(i).get()->second == (i+3.0));
    }

    world.gop.fence();
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        test0(world);
        test1(world);
        test1(world);
        test1(world);
    }
    catch (MPI::Exception e) {
        error("caught an MPI exception");
    }
    catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}
