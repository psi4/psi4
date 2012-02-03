#ifndef _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_
#define _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_

#include <psi4-dec.h>
#include <map>

namespace psi{

  /**
   * The MOSpace class is used to define orbital spaces in which to transform
   * integrals
   */
class MOSpace{

    public:
        ~MOSpace();
        MOSpace(const char label, const int nirreps, const int *aOrbsPI,
                const int *bOrbsPI, const double ***aEvecs, const double ***bEvecs,
                const int *aIndex, const int *bIndex);
        /**
         * The MOSpace::frc space can be used to define the frozen occupied space.
         *
         * The label associated with this space is 'o'
         */
        #define MOSPACE_FZC 'o'
        static boost::shared_ptr<MOSpace> fzc;
        /**
         * The MOSpace::occ space can be used to define the occupied space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  For restricted transformations, this corresponds to
         * singly- plus doubly-occupied orbitals, for unrestricted transforms only
         * occupied orbitals are included
         *
         * The label associated with this space is 'O'
         */
        #define MOSPACE_OCC 'O'
        static boost::shared_ptr<MOSpace> occ;
        /**
         * The MOSpace::frv space can be used to define the frozen virtual space.
         *
         * The label associated with this space is 'v'
         */
        #define MOSPACE_FZV 'v'
        static boost::shared_ptr<MOSpace> fzv;
        /**
         * The MOSpace::vir space can be used to define the virtual space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  For restricted transformations, this corresponds to
         * singly occupied plus virtual orbitals, for unrestricted transforms only
         * virtual orbitals are included
         *
         * The label associated with this space is 'V'
         */
        #define MOSPACE_VIR 'V'
        static boost::shared_ptr<MOSpace> vir;
        /**
         * The MOSpace::all space can be used to define the full MO space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  All active molecular orbtitals are transformed
         *
         * The label associated with this space is 'a'
         */
        #define MOSPACE_ALL 'A'
        static boost::shared_ptr<MOSpace> all;
        /**
         * The MOSpace::nil space can be used to define the atomic orbital space.
         *
         * The label associated with this space is 'n'
         */
        #define MOSPACE_NIL 'n'
        static boost::shared_ptr<MOSpace> nil;

        // These are to allow the map to be used
        friend bool operator==(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label == rhs._label; }
        friend bool operator>=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label >= rhs._label; }
        friend bool operator!=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label != rhs._label; }
        friend bool operator<=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label <= rhs._label; }
        friend bool operator<(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label < rhs._label; }
        friend bool operator>(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs._label > rhs._label; }
        friend bool operator==(const MOSpace &lhs, const char c)
                                { return lhs._label == c; }
        friend bool operator>=(const MOSpace &lhs, const char c)
                                { return lhs._label >= c; }
        friend bool operator!=(const MOSpace &lhs, const char c)
                                { return lhs._label != c; }
        friend bool operator<=(const MOSpace &lhs, const char c)
                                { return lhs._label <= c; }
        friend bool operator<(const MOSpace &lhs, const char c)
                                { return lhs._label < c; }
        friend bool operator>(const MOSpace &lhs, const char c)
                                { return lhs._label > c; }

        /// Get the unique identifier for this space
        char label() {return _label;}

        /// Get the number of alpha orbitals per irrep
        const int* aOrbsPI() {return _aOrbsPI;}

        /// Get the number of beta orbitals per irrep
        const int* bOrbsPI() {return _bOrbsPI == NULL ? _aOrbsPI : _bOrbsPI;}

        /// Get the alpha orbital symmetries
        const int* aOrbSym() {return _aOrbSym;}

        /// Get the beta orbital symmetries
        const int* bOrbSym() {return _bOrbSym == NULL ? _aOrbSym : _bOrbSym;}

        /// Get the alpha orbital indexing array for IWL
        const int* aIndex() {return _aIndex;}

        /// Get the beta orbital indexing array for IWL
        const int* bIndex() {return _bIndex == NULL ? _aIndex : _bIndex;}

        /// Get the alpha MO coefficients for this space
        const double*** aEvecs() {return _aEvecs;}

        /// Get the beta MO coefficients for this space
        const double*** bEvecs() {return _bEvecs == NULL ? _aEvecs : _bEvecs;}

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
        // The number of alpha orbitals per irrep
        const int *_aOrbsPI;
        // The number of beta orbitals per irrep
        const int *_bOrbsPI;
        // The alpha indexing array
        const int *_aIndex;
        // The beta indexing array
        const int *_bIndex;
        // The symmetry of each alpha orbital
        int *_aOrbSym;
        // The symmetry of each beta orbital
        int *_bOrbSym;
        // The alpha eigenvector for this space
        const double ***_aEvecs;
        // The beta eigenvector for this space
        const double ***_bEvecs;
};

} // End namespaces

#endif // Header guard
