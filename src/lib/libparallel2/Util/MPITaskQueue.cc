/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include "MPITaskQueue.h"
namespace psi{
namespace LibParallel{
int MPITaskQueue::operator[](const int i)const{
   return Tasks_[i];
}

MPITaskQueue::MPITaskQueue():LastTask_(0){

}

void MPITaskQueue::Copy(const MPITaskQueue& other){
   this->Tasks_=other.Tasks_;
   this->LastTask_=other.LastTask_;
   this->Scheduler_=other.Scheduler_;
}

int MPITaskQueue::size()const{
   return Tasks_.size();
}

MPITaskQueue::MPITaskQueue(const MPITaskQueue& other){
   this->Copy(other);
}

const MPITaskQueue& MPITaskQueue::operator=(const MPITaskQueue& other){
   if(this!=&other)this->Copy(other);
   return *this;
}

void MPITaskQueue::operator+=(const MPITaskQueue& other){
   this->Tasks_.reserve(Tasks_.size()+other.Tasks_.size());
   this->Tasks_.insert(Tasks_.end(),other.Tasks_.begin(),other.Tasks_.end());
}

MPITaskQueue& MPITaskQueue::operator<<(const int NewTask){
   this->Tasks_.push_back(NewTask);
   //Make sure our new task will actually be returned when the next
   //next is called
   if(LastTask_==this->size()-1)LastTask_--;
   return *this;
}

MPITaskQueue MPITaskQueue::operator+(const MPITaskQueue& other){
   MPITaskQueue temp(*this);
   temp+=other;
   return temp;
}

int MPITaskQueue::Begin(){
   LastTask_=0;
   return Tasks_[LastTask_];
}

int MPITaskQueue::Next(){
   ++LastTask_;
   //For the last iteration can't increment LastTask_ or we segfault
   int next=(LastTask_==Tasks_.size()? -1 : LastTask_);
   return (next==-1?next:Tasks_[next]);
}

bool MPITaskQueue::Done()const{
   return (LastTask_==Tasks_.size());
}
}}//End namespaces