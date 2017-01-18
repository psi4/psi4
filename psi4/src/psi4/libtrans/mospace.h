/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_
#define _PSI_SRC_LIB_LIBTRANS_MOSPACE_H_

#include "psi4/psi4-dec.h"
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
        MOSpace(const char label, const std::vector<int> orbsPI);

        /**
         * The MOSpace::frc space can be used to define the frozen occupied space.
         *
         * The label associated with this space is 'o'
         */
        #define MOSPACE_FZC 'o'
        static std::shared_ptr<MOSpace> fzc;
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
        static std::shared_ptr<MOSpace> occ;
        /**
         * The MOSpace::frv space can be used to define the frozen virtual space.
         *
         * The label associated with this space is 'v'
         */
        #define MOSPACE_FZV 'v'
        static std::shared_ptr<MOSpace> fzv;
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
        static std::shared_ptr<MOSpace> vir;
        /**
         * The MOSpace::all space can be used to define the full MO space.  Frozen
         * orbitals are handled consistently with how the transformation object is
         * constructed.  All active molecular orbtitals are transformed
         *
         * The label associated with this space is 'A'
         */
        #define MOSPACE_ALL 'A'
        static std::shared_ptr<MOSpace> all;
        /**
         * The MOSpace::nil space can be used to define the atomic orbital space.
         *
         * The label associated with this space is 'n'
         */
        #define MOSPACE_NIL 'n'
        static std::shared_ptr<MOSpace> nil;
        /**
         * The MOSpace::dum space is a dummy space with a single function in each irrep.
         * It is used for converting a single aux index into a DPD pair.
         *
         * The label associated with this space is 'd'
         */
        #define MOSPACE_DUM 'd'
        static std::shared_ptr<MOSpace> dum;

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
        const std::vector<int>& bIndex() {return bIndex_.size() == 0 ? aIndex_ : bIndex_;}

        /// Whether this is just a placeholder
        bool placeholder() const {return placeholder_;}
        void set_placeholder(bool t_f) {placeholder_ = t_f;}


    protected:
        MOSpace(const char label);
        /**
         * The identifier for this space; this must be unique for each space.  See the
         * documentation for the static spaces in MOSpace to see which labels have already
         * been used
         */
        const char label_;
        // The indices (Pitzer) of the alpha orbitals
        std::vector<int> aOrbs_;
        // The indices (Pitzer) of the beta orbitals
        std::vector<int> bOrbs_;
        // The alpha reindexing array
        std::vector<int> aIndex_;
        // The beta reindexing array
        std::vector<int> bIndex_;
        // Whether this only describes dimensions, and has no orbitals associated with it
        bool placeholder_;
};

} // End namespaces

#endif // Header guard
