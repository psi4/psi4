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
      FragmentedSystem(const SuperCell& System2Frag, const int N=1);
      std::string PrintOut(const int Value=1)const;
      typedef std::vector<boost::shared_ptr<Fragment> >::const_iterator
            iterator;
      int N()const;
      double Coef(const int N, const int i)const;
      iterator begin(const int N)const;
      iterator end(const int N)const;
};


}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_FRAGMENTEDSYS_H_ */
