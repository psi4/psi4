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

#include "mospace.h"
#include "psi4/libciomr/libciomr.h"
#include <vector>

using namespace psi;

/// Keeps track of which labels have been assigned, for safety
std::shared_ptr<MOSpace> MOSpace::fzc(new MOSpace(MOSPACE_FZC));
std::shared_ptr<MOSpace> MOSpace::occ(new MOSpace(MOSPACE_OCC));
std::shared_ptr<MOSpace> MOSpace::vir(new MOSpace(MOSPACE_VIR));
std::shared_ptr<MOSpace> MOSpace::fzv(new MOSpace(MOSPACE_FZV));
std::shared_ptr<MOSpace> MOSpace::all(new MOSpace(MOSPACE_ALL));
std::shared_ptr<MOSpace> MOSpace::nil(new MOSpace(MOSPACE_NIL));
std::shared_ptr<MOSpace> MOSpace::dum(new MOSpace(MOSPACE_DUM));

/**
 * This creates an empty MOSpace with just a label.  This is solely for the
 * construction of the pre-defined spaces; use the longer form of the constructor
 * for custom spaces.
 */
MOSpace::MOSpace(char label):
        label_(label),
        aOrbs_(0),
        bOrbs_(0),
        aIndex_(0),
        bIndex_(0)
{
}

/**
 * Defines a custom orbital space with different alpha and beta spaces
 * @param label  - a single character to label this space.  This must be unique to this
 *                 space, so see the MOSpace static member variables for a list of the
 *                 labels already used.  The uniqueness is checked internally.
 * @param aOrbs  - an array of dimension <= nso, describing the Pitzer indices of the alpha
 *                 (and beta) orbitals present
 * @param aIndex - an array of dimension #orbitals describing the number of each alpha (and beta)
 *                 orbital in the space. This is only for the purposes of IWL output, so it can
 *                 be passed as an empty vector for DPD output.
 */
MOSpace::MOSpace(const char label, const std::vector<int> aOrbs, const std::vector<int> aIndex):
        label_(label),
        aOrbs_(aOrbs),
        bOrbs_(aOrbs),
        aIndex_(aIndex),
        bIndex_(aIndex),
        placeholder_(false)
{
}

/**
 * Defines a custom orbital space with different alpha and beta spaces
 * @param label  - a single character to label this space.  This must be unique to this
 *                 space, so see the MOSpace static member variables for a list of the
 *                 labels already used.  The uniqueness is checked internally.
 * @param aOrbs  - an array of dimension <= nso, describing the Pitzer indices of the alpha orbitals present
 * @param bOrbs  - an array of dimension <= nso, describing the Pitzer indices of the beta orbitals present
 * @param aIndex - an array of dimension #orbitals describing the number of each alpha orbital in
 *                 the space. This is only for the purposes of IWL output, so it can be passed as
 *                 an empty vector for DPD output.
 * @param bIndex - an array of dimension #orbitals describing the number of each beta orbital in
 *                 the space. For restricted transformations or for DPD output only, this can be
 *                 passed as an empty vector.
 */
MOSpace::MOSpace(const char label, const std::vector<int> aOrbs, const std::vector<int> bOrbs,
                 const std::vector<int> aIndex, const std::vector<int> bIndex):
        label_(label),
        aOrbs_(aOrbs),
        bOrbs_(bOrbs),
        aIndex_(aIndex),
        bIndex_(bIndex),
        placeholder_(false)
{
}

/**
 * Defines a custom orbital space over (symmetrized) atomic orbitals.  Used for auxilliary basis sets.
 * @param label  - a single character to label this space.  This must be unique to this
 *                 space, so see the MOSpace static member variables for a list of the
 *                 labels already used.  The uniqueness is checked internally.
 * @param orbsPI - the numer of orbitals per irrep in this SO space
 */
MOSpace::MOSpace(const char label, const std::vector<int> orbsPI):
        label_(label),
        aOrbs_(orbsPI),
        placeholder_(true)
{
}


MOSpace::~MOSpace()
{
}
