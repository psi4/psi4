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

#include "MoleculeGuts.h"
#include "Implementations/MolItrGuts.h"
#include "../Atom.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<const Atom> SharedAtom;
typedef boost::shared_ptr<MolItrGuts> SharedItr;
MolItr MoleculeGuts::Begin()const{
   return MolItr(SharedItr(new MolItrGuts(Atoms_.begin())));
}

MolItr MoleculeGuts::End()const{
   return MolItr(SharedItr(new MolItrGuts(Atoms_.end())));
}
void MoleculeGuts::Copy(const MoleculeGuts& other){
   this->Mult_=other.Mult_;
   this->Charge_=other.Charge_;
   MolItr Start=other.Begin(),End=other.End();
   std::vector<SharedAtom> Temp;
   for(;Start!=End;++Start)
      Temp.push_back(SharedAtom(new Atom((*(*Start)))));
   this->Atoms_=Temp;

}

}}//End namespaces

