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

// User defined key class that we don't want to modify
struct Key {
    int k;
    Key() : k(-1) {}
    Key(int k) : k(k) {}

    bool operator==(const Key& b) const {
        return k==b.k;
    }
};

ostream& operator<<(ostream&s, const Key& key) {
    s << "Key(" << key.k << ")";
    return s;
}

// Make the key serialiable using non-intrusive mechanism
namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveSerializeImpl<Archive,Key> {
            static inline void serialize(const Archive& ar, Key& obj) {
                ar & obj.k;
            }
        };
    }
}

// Make the key hashable using non-intrusive mechanism
hashT hash_value(const Key& key) {
    return key.k;
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    WorldContainer<Key,double> fred(world);

    fred.replace(Key(99),99.0);

    cout << fred.find(Key(99)).get()->second << endl;

    WorldContainer<Key,double>::iterator it = fred.find(Key(99));
    // WorldContainer<Key,double>::pairT& p = *it;
    // cout << p;
    cout << *it;

    WorldContainer<Key,double>::const_iterator c_it = it;
    const WorldContainer<Key,double>::pairT& cp = *c_it;
    cout << cp;

    // This fails because as shown above cannot derefence const iterator
    // ... works OK for non-const
    cout << *it << endl;
    cout << *c_it << endl;

    finalize();
    return 0;
}
