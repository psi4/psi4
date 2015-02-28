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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENTEDSYS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENTEDSYS_H_

#include "Implementations/FragSysGuts.h"
namespace psi{
namespace LibMolecule{

class FragmentedSystem:private FragSysGuts{
   public:
      ///Fragments a molecule
      FragmentedSystem(const Molecule& System2Frag,const int N=1);
      ///Special constructor for supercells,
      FragmentedSystem(const SuperCell& System2Frag, const int N=1);
      ///Prints out the fragments in a pretty format
      std::string PrintOut(const int Value=1)const;
      ///A typedef of an iterator over fragments
      typedef std::vector<boost::shared_ptr<Fragment> >::const_iterator
            iterator;
      ///Returns the overall MBE order
      int N()const;
      ///Returns the coefficient of N-mer i
      double Coef(const int N, const SerialNumber& i)const;
      ///Returns an iterator to the first N-Mer
      iterator begin(const int N)const;
      ///Returns an iterator just past the last N-Mer
      iterator end(const int N)const;
      /** \brief Returns the serial number of the n-mer that has the property
       *         of interest.
       *
       *  Consider a symmetric system in which dimer2 is a pure rotation plus
       *  translation of dimer1.  In this case we only want to compute that
       *  property once, say for dimer1, and then multiply it by 2.
       *  In this scenario this function would return the serial number of
       *  dimer1, when it is given dimer2.  When given dimer1 it simply
       *  returns the serial number of dimer1.  The returned value
       *  will be how your property is indexed.
       */
      const SerialNumber& SNLookUp(const SerialNumber& SN)const;
};
}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENTEDSYS_H_ */
