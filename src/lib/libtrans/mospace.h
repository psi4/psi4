#ifndef _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_
#define _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_

#include <psi4-dec.h>
#include <map>

using namespace psi;

namespace psi{ namespace libtrans{

class SpaceInfo;

  /**
   * The MOSpace class is used to define orbital spaces in which to transform
   * integrals
   */
class MOSpace{

    public:
        ~MOSpace();
        MOSpace(const char label, int *aFOrbPI, int *bFOrbPI, int *aOrbsPI, int *bOrbsPI);
        /**
         * The MOSpace::occ space can be used to define the occupied space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  For restricted transformations, this corresponds to
         * singly- plus doubly-occupied orbitals, for unrestricted transforms only
         * occupied orbitals are included
         *
         * The label associated with this space is 'o'
         */
        #define MOSPACE_OCC 'o'
        static shared_ptr<MOSpace> occ;
        /**
         * The MOSpace::vir space can be used to define the virtual space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  For restricted transformations, this corresponds to
         * singly occupied plus virtual orbitals, for unrestricted transforms only
         * virtual orbitals are included
         *
         * The label associated with this space is 'v'
         */
        #define MOSPACE_VIR 'v'
        static shared_ptr<MOSpace> vir;
        /**
         * The MOSpace::all space can be used to define the virtual space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  All active molecular orbtitals are transformed
         *
         * The label associated with this space is 'a'
         */
        #define MOSPACE_ALL 'a'
        static shared_ptr<MOSpace> all;
        /**
         * The MOSpace::nil space can be used to define the atomic orbital space.
         *
         * The label associated with this space is 'n'
         */
        #define MOSPACE_NIL 'n'
        static shared_ptr<MOSpace> nil;
        
        // Returns the unique identifier for this space
        const char label() {return _label;}

        // These are to allow the map to be used
        friend bool operator==(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() == rhs.label(); }
        friend bool operator>=(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() >= rhs.label(); }
        friend bool operator!=(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() != rhs.label(); }
        friend bool operator<=(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() <= rhs.label(); }
        friend bool operator<(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() < rhs.label(); }
        friend bool operator>(MOSpace &lhs, MOSpace &rhs)
                                { return lhs.label() > rhs.label(); }
        friend bool operator==(MOSpace &lhs, char c)
                                { return lhs.label() == c; }
        friend bool operator>=(MOSpace &lhs, char c)
                                { return lhs.label() >= c; }
        friend bool operator!=(MOSpace &lhs, char c)
                                { return lhs.label() != c; }
        friend bool operator<=(MOSpace &lhs, char c)
                                { return lhs.label() <= c; }
        friend bool operator<(MOSpace &lhs, char c)
                                { return lhs.label() < c; }
        friend bool operator>(MOSpace &lhs, char c)
                                { return lhs.label() > c; }
    protected:
        MOSpace(const char label);
        /**
         * The identifier for this space; this must be unique for each space.  See the
         * documentation for the static spaces in MOSpace to see which labels have already
         * been used
         */
        const char _label;
        // This keeps track of which labels have been assigned by other spaces
        static std::map<char, int> labelsUsed;
        // The indexing and orbitals for this spaces
        shared_ptr<SpaceInfo> _spaceInfo;
};

}} // End namespaces

// This is here so that files including this have clean(er) syntax
using namespace psi::libtrans;

#endif // Header guard