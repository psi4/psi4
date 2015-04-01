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

#include <sstream>
#include "LibFragMolecule.h"
#include "Atom.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<Atom> SharedAtom;
typedef std::vector<SharedAtom>::iterator AtomItr;
MolItr Molecule::Begin()const{
   return MoleculeGuts::Begin();
}

MolItr Molecule::End()const{
   return MoleculeGuts::End();
}

const Molecule& Molecule::operator+=(const Molecule& other){
   MolItr Start=other.Begin(),Done=other.End();
   for(;Start!=Done;++Start)(*this)<<(*(*Start));
   return *this;
}

Molecule Molecule::operator+(const Molecule& other){
   Molecule temp(*this);
   temp+=other;
   return temp;
}


std::string Molecule::PrintOut(const int DebugLevel)const{
   std::stringstream message;
   MolItr Atoms=Begin(),AtomsEnd=End();
   for(;Atoms!=AtomsEnd;++Atoms)
      message<<Atoms->PrintOut(false,DebugLevel)
             <<std::endl;
   return message.str();
}

double Molecule::operator()(const int AtomI,const int i)const{
   return (*LookUp(AtomI))[i];
}

bool Molecule::operator==(const Molecule& other)const{
   if(this->NAtoms()!=other.NAtoms())return false;
   MolItr AtomI=this->Begin(),AtomIEnd=this->End();
   for(int i=0;AtomI!=AtomIEnd;++AtomI,++i){
      MolItr AtomJ=other.Begin(),AtomJEnd=other.End();
      bool AreEqual=false;
      for(;AtomJ!=AtomJEnd;++AtomJ)
         if((*(*AtomI))==(*(*AtomJ))){
            AreEqual=true;
            break;
         }
      if(!AreEqual)return false;
   }
   return true;
}

}}//end namespaces

