/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */
#ifndef PSIAPI_FRAGMENT_H_
#define PSIAPI_FRAGMENT_H_

namespace PsiAPI{

/** \brief A subset (not necessarily "proper") of a Molecule.
 *
 *   don't want an atom to always be a ghost and we don't want
 *   to have a ton of copies sitting around (particularly bad for the MBE).
 *   In these cases we have the following filling paradigm, which
*   for concreteness we limit ourselves to a dimer CP correction, but the
*   generalization should be obvious.  Assume we have some dimer \f$AB\f$,
*   first we create a molecule \f$M\f$ which is to serve as the "universe"
*   in a set sense.  We define three subsets of \f$M\f$: \f$AB\f$, \f$A\f$
*   and \f$B\f$, which we call fragments in keeping with the traditional
*   nomenclature of electronic structure theory.
*
*
*   The general strategy is: we fill \f$M\f$ with the atoms in \f$A\f$
*   while simultaneously adding them to \f$A\f$, repeat for \f$B\f$,
*   fill \f$AB\f$ with the union of \f$A\f$ and \f$B\f$, add copies of
*   \f$B\f$ as ghosts to \f$M\f$ and \f$A\f$, finally repeat for \f$A\f$
*   as ghosts to \f$B\f$.  In doing this
*
*
*   Code-wise we have:
*
    \code
    //Our "universe"
    SharedPtr<PsiAPI::Molecule> M;
    //Our fragments
    PsiAPI::Fragment AB(M),A(M),B(M);
    //The same molecule from the last code example (presumed filled)
    PsiAPI::Molecule MyMol;
    //A set of pointers to the atoms in MyMol that are in A
    std::set<SharedPtr<Atom> > InA;


    \endcode
class Fragment{

};

}//End namespace


#endif /* PSIAPI_FRAGMENT_H_ */
