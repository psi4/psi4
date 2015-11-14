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
#include <iostream>
#include "Molecule.h"
#include "Atom.h"
namespace PsiAPI{

void Molecule::AddAtom(size_t Z,double x,double y, double z, double Q, double M){
   //Need to watch for the vector reallocating memory, which invalidates
   //references
   bool reallocate=false;
   //Reallocation is going to occur on one of the push_back's
   if(Carts_.capacity()-Carts_.size()<3)reallocate=true;
   Carts_.push_back(x);
   Carts_.push_back(y);
   Carts_.push_back(z);
   double* Start=&(Carts_[Carts_.size()-3]);
   if(reallocate){
      std::vector<SharedAtom_t> Temp;
      Temp.reserve(Carts_.size()/3);
      cAtomItr AtomI=this->AtomBegin(),AtomIEnd=this->AtomEnd();
      for(size_t counter=0;AtomI!=AtomIEnd;++AtomI,counter+=3){
         Atom AI=(*(*AtomI));
         Temp.push_back(SharedAtom_t(
               new Atom(AI.Z(),&Carts_[counter],AI.Q(),AI.Mass())));
      }
      Atoms_=Temp;
   }
   Atoms_.push_back(SharedAtom_t(new Atom(Z,Start,Q,M)));
}


Molecule::Molecule(double Charge,size_t Mult,size_t NAtoms):
   Charge_(Charge),Mult_(Mult){
   if(NAtoms>0){
      Carts_.reserve(3*NAtoms);
      Atoms_.reserve(NAtoms);
   }
}


Molecule::operator std::string()const{
   cAtomItr AtomI=AtomBegin(),AtomIEnd=AtomEnd();
   std::stringstream Result;
   for(;AtomI!=AtomIEnd;++AtomI)
      Result<<std::string(*(*AtomI))<<std::endl;
   return Result.str();
}

}


