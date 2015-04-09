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

#include "MolItrGuts.h"

namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<IteratorGuts> SharedGuts;
SharedGuts MolItrGuts::Clone()const{
   SharedGuts temp(new MolItrGuts(*this));
   return temp;
}

MolItrGuts::MolItrGuts(const MolItrGuts& other):IteratorGuts(other){
   this->Copy(other);
}


bool MolItrGuts::IsEqual(const IteratorGuts& other)const{
   return this->MyItr_==UpCast(other)->MyItr_;
}

MolItrGuts::MolItrGuts(const Itr_t& MyItr):MyItr_(MyItr){}

void MolItrGuts::Copy(const MolItrGuts& other){
   this->MyItr_=other.MyItr_;
}


const MolItrGuts& MolItrGuts::operator=(const MolItrGuts& other){
   if(this!=&other)this->Copy(other);
   return *this;
}

const MolItrGuts* MolItrGuts::UpCast(const IteratorGuts& other)const{
   return dynamic_cast<const MolItrGuts*>(& other);
}

}}//End namespaces

