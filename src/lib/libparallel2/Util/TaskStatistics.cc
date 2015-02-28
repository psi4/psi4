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

#include "TaskStatistics.h"
#include "../Communicator.h"
#include "../TaskJobGuts/MPITaskGuts.h"
namespace psi{
namespace LibParallel{
typedef boost::shared_ptr<std::vector<MPITaskGuts> > SharedTasks;

int TaskStats::NPriorities()const{return PriBins_.size();}

int TaskStats::NTimes(const int i)const{return PriBins_[i];}

SchedAlgorithm TaskStats::SuggestedAlgorithm()const{return Sugg_;}

TaskStats::TaskStats(SharedTasks& Tasks,const bool ForceDyn,
      boost::shared_ptr<const Communicator>& State):
      Sugg_(SIMPLE),UnknownPri_(false){
    for(int i=0;i<Tasks->size();i++){
       int Prior=(*Tasks)[i].Priority();
       if(Prior<0){
          UnknownPri_=true;
          Prior=-1;
       }
       if(PriBins_.count(Prior)==0)PriBins_[Prior]=1;
       else PriBins_[Prior]++;
    }
    if(UnknownPri_)HandleUnknown(Tasks,State);
    if(ForceDyn)Sugg_=DYNAMICRR;
    else if(NPriorities()==1&&!UnknownPri_)Sugg_=SIMPLE;
    else Sugg_=STATICRR;
}

void TaskStats::HandleUnknown(
      boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks,
      boost::shared_ptr<const Communicator>& State){
   int NUnknown=this->NTimes(-1);
   std::vector<int> NewPrior(Tasks->size());
   for(int i=0;i<Tasks->size();i++)NewPrior[i]=this->RdmNum(10);
   State->Bcast(&NewPrior[0],Tasks->size(),0);
   for(int i=0;i<Tasks->size();i++)(*Tasks)[i].SetPriority(NewPrior[i]);
   if(State->Me()==0){
      (*this)<<"Printing in Handle Unknown!!!\n";
      for(int i=0;i<Tasks->size();i++)(*Tasks)[i].PrintOut();
   }
}

}}

