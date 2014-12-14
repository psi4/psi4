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

#include "MPITaskGuts.h"

#include <iostream>
namespace psi{
namespace LibParallel{

void MPITaskGuts::SetPriority(const int i){
   Task_.first=i;
}

void MPITaskGuts::PrintOut()const{
   (*this)<<"Priority= "<<this->Priority();
   (*this)<<" OrigNum= "<<this->Number()<<"\n";
}


void MPITaskGuts::Copy(const MPITaskGuts& other){
   this->Task_=other.Task_;
}

MPITaskGuts::MPITaskGuts(const MPITaskGuts& other){
   this->Copy(other);
}

const MPITaskGuts& MPITaskGuts::operator=(const MPITaskGuts& other){
   if(this!=&other)this->Copy(other);
   return *this;
}

MPITaskGuts::MPITaskGuts(const int Priority,const int Number){
   Task_.first=Priority;
   Task_.second=Number;
}

int MPITaskGuts::Priority()const{return Task_.first;}


int MPITaskGuts::Number()const{return Task_.second;}

bool MPITaskGuts::operator<(const MPITaskGuts& other)const{
   return (this->Task_.first<other.Task_.first);
}

bool MPITaskGuts::operator==(const MPITaskGuts& other)const{
   return (this->Task_.first==other.Task_.first);
}

bool MPITaskGuts::operator<=(const MPITaskGuts& other)const{
   return ((*this)<other?true:(*this)==other);
}

bool MPITaskGuts::operator>(const MPITaskGuts& other)const{
   return (!((*this)<=other));
}

bool MPITaskGuts::operator>=(const MPITaskGuts& other)const{
   return (!((*this)<other));
}

}}


