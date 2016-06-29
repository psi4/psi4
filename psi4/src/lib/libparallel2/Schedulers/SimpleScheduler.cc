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

#include "SimpleScheduler.h"
#include "../Util/MPITaskQueue.h"
namespace psi{
namespace LibParallel{
SimpleScheduler::SimpleScheduler(const int NTasks,
      boost::shared_ptr<const Communicator> State) :
      StaticScheduler(State),
      NTasks_(NTasks), Remainder_(0),
      BatchSize_(0), MyStart_(0), MyEnd_(0){
      Alg_=SIMPLE;
      SetValues();
}


void SimpleScheduler::SetValues(){
   int NProc=IState_->NProc();
   int Me=IState_->Me();
   Remainder_=NTasks_%NProc;
   BatchSize_=(NTasks_-Remainder_)/NProc;
   MyStart_=Me*BatchSize_;
   MyEnd_=(Me+1)*BatchSize_;
   for(int i=MyStart_;i<MyEnd_;i++)(*Tasks_)<<i;
   if(Me<Remainder_)(*Tasks_)<<(BatchSize_*NProc+Me);
}

}}//End namespaces