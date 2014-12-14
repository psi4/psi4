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

#include "Communicator.h"
#include "ParallelEnvironment.h"
#include <iostream>
namespace psi{
namespace LibParallel{
typedef boost::shared_ptr<Communicator> SharedThis;

void Communicator::PrintOut()const{
   (*this)<<"I am process "<<Me()<<"/"<<NProc()<<
           " on comm: "<<Name()<<"\n";
}

bool Communicator::Active()const{return CommGuts::Active();}

Communicator::Communicator(const Communicator& other):CommGuts(other){}
const Communicator& Communicator::operator=(const Communicator& other){
   CommGuts::operator=(other);
   return *this;
}
Communicator::Communicator(const std::string& Name,
      ParallelEnvironmentGuts* Parent):CommGuts(Name,Parent){}

SharedThis Communicator::MakeComm(const int Color,
      const std::string& Name)const{
   if(!this->Active())this->Error("Current communicator has been freed.  "
         "You can't use it to make another communicator");
   boost::shared_ptr<CommGuts> temp=CommGuts::MakeComm(Name,Color);
   SharedThis temp2(new Communicator(*this));
   temp2->CommGuts::operator=((*temp));
   CommGuts::RegisterComm(temp2);
   return temp2;
}
int Communicator::Me()const{return CommGuts::Me();}
int Communicator::NProc()const{return CommGuts::NProc();}
void Communicator::FreeComm(){return CommGuts::FreeComm();}
void Communicator::Barrier()const{CommGuts::Barrier();}
int Communicator::Probe(const int Sender,const int MessageTag,
                        bool Block)const{
   return CommGuts::Probe(Sender,MessageTag,Block);
}

}}//End namespaces

