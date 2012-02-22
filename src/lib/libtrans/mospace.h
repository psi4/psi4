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
        MOSpace(const char label, const std::vector<int> aOrbs, const std::vector<int> bOrbs,
                const std::vector<int> aIndex, const std::vector<int> bIndex);
        MOSpace(const char label, const std::vector<int> aOrbs, const std::vector<int> aIndex);

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
                                { return lhs.label_ == rhs.label_; }
        friend bool operator>=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs.label_ >= rhs.label_; }
        friend bool operator!=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs.label_ != rhs.label_; }
        friend bool operator<=(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs.label_ <= rhs.label_; }
        friend bool operator<(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs.label_ < rhs.label_; }
        friend bool operator>(const MOSpace &lhs, const MOSpace &rhs)
                                { return lhs.label_ > rhs.label_; }
        friend bool operator==(const MOSpace &lhs, const char c)
                                { return lhs.label_ == c; }
        friend bool operator>=(const MOSpace &lhs, const char c)
                                { return lhs.label_ >= c; }
        friend bool operator!=(const MOSpace &lhs, const char c)
                                { return lhs.label_ != c; }
        friend bool operator<=(const MOSpace &lhs, const char c)
                                { return lhs.label_ <= c; }
        friend bool operator<(const MOSpace &lhs, const char c)
                                { return lhs.label_ < c; }
        friend bool operator>(const MOSpace &lhs, const char c)
                                { return lhs.label_ > c; }

        /// Get the unique identifier for this space
        char label() {return label_;}

        /// Get the alpha orbitals
        const std::vector<int>& aOrbs() const {return aOrbs_;}

        /// Get the beta orbitals
        const std::vector<int>& bOrbs() const {return bOrbs_;}

        /// Get the alpha orbital indexing array for IWL
        const std::vector<int>& aIndex() {return aIndex_;}

        /// Get the beta orbital indexing array for IWL
        const std::vector<int>& bIndex() {return bIndex_ == NULL ? aIndex_ : bIndex_;}


    protected:
        MOSpace(const char label);
        /**
         * The identifier for this space; this must be unique for each space.  See the
         * documentation for the static spaces in MOSpace to see which labels have already
         * been used
         */
        const char label_;
        // This keeps track of which labels have been assigned by other spaces
        static std::map<char, int> labelsUsed;
        // The indices (Pitzer) of the alpha orbitals
        std::vector<int> aOrbs_;
        // The indices (Pitzer) of the beta orbitals
        std::vector<int> bOrbs_;
        // The alpha reindexing array
        std::vector<int> aIndex_;
        // The beta reindexing array
        std::vector<int> bIndex_;
};

} // End namespaces

#endif // Header guard
