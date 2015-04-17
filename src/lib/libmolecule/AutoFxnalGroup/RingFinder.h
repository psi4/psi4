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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_
#include "LinearSearch.h"

namespace psi{
namespace LibMolecule{

/** \brief The class that finds rings
 *
 *   This class is basically a linear search except that it ensures that
 *   the groups on the ends of the ring are attached.
 *
 */
template<typename...Groups>
class RingFinder: public LinearSearch<Groups...>{
   private:
      typedef LinearSearch<Groups...> Base_t;
   protected:
      std::vector<size_t> DeterminePrior(
            std::deque<boost::shared_ptr<Node> >&,
            bool& Symm)const;
   public:
      RingFinder(const std::string& BaseAbbrv,const std::string& BaseName):
                    Base_t(BaseAbbrv,BaseName){}
      bool FindMe(boost::shared_ptr<Node> Nodes);
};


/*****************Implementations************************/

template<typename...Groups>
bool RingFinder<Groups...>::FindMe(boost::shared_ptr<Node> Nodes){
   typedef boost::shared_ptr<Node> SharedNode;
   //Can't possibly make a ring
   if(Nodes->NEdges()<2)return false;
   if(!Base_t::FindMe(Nodes))return false;
   //Check if NodeI is connected to NodeJ
   SharedNode NodeI,NodeJ;
   Node::iterator It=this->SubBegin(),ItEnd=this->SubEnd();
   NodeI=*It;
   for(;It!=ItEnd;++It)NodeJ=*It;
   Node::ConnItr It1=NodeI->ConnBegin(),It1End=NodeI->ConnEnd();
   for(;It1!=It1End;++It1)if(It1->second.get()==NodeJ.get())return true;
   return false;
}

template<typename...Groups>
std::vector<size_t> RingFinder<Groups...>::DeterminePrior(
      std::deque<boost::shared_ptr<Node> >& FN,
      bool& Symm)const{
   typedef boost::shared_ptr<Node> SharedNode;
   typedef std::deque<SharedNode> FN_t;
   typedef FN_t::iterator cFN_Itr;
   typedef std::vector<size_t> Pri_t;
   typedef Node::const_iterator NodeItr;
   cFN_Itr It=FN.begin(),ItEnd=FN.end();
   Pri_t Priors,Subs;
   FN_t Order=FN;
   for(;It!=ItEnd;++It){
      FN_t TrialLayout;
      TrialLayout.push_back(*It);
      cFN_Itr It1=It;++It1;
      for(;It1!=ItEnd;++It1)TrialLayout.push_back(*It1);
      for(It1=FN.begin();It1!=It;++It1)TrialLayout.push_back(*It1);
      Pri_t TempPrior;
      TempPrior=Node::DeterminePrior(TrialLayout);
      It1=TrialLayout.begin();
      Pri_t TempSubs;
      cFN_Itr It1End=TrialLayout.end();
      for(size_t counter=0;It1!=It1End;++It1,++counter){
         Node::ConnItr It2=(*It1)->ConnBegin(),It2End=(*It1)->ConnEnd();
         for(;It2!=It2End;++It2){
           cFN_Itr It3=TrialLayout.begin();
           for(;It3!=It1End;++It3)if(It3->get()==It2->second.get())break;
           if(It3!=It1End)continue;
           TempSubs.push_back(counter);
         }
      }
      if(TempSubs.size()==0){
         Subs=TempSubs;
         Order=TrialLayout;
         Priors=TempPrior;
         break;
      }
      if(TempSubs[0]!=0)continue;
      bool good=false;

      if(Subs.size()!=0){
         Pri_t::iterator Pi=TempSubs.begin(),PiEnd=TempSubs.end(),
               Pj=Subs.begin();
         for(;Pi!=PiEnd;++Pi,++Pj){
            if(*Pj<*Pi)break;
            if(*Pi<*Pj){
               good=true;
               break;
            }
         }
         if(Pi!=PiEnd&&!good)continue;
         if(!good){
            Pi=TempPrior.begin();PiEnd=TempPrior.end();Pj=Priors.begin();
            for(;Pi!=PiEnd;++Pi,++Pj){
               if(*Pj<*Pi)break;
               if(*Pi<*Pj){
                  good=true;
                  break;
               }
            }
         }
         if(!good)continue;
      }
      Subs=TempSubs;
      Order=TrialLayout;
      Priors=TempPrior;
   }
   FN=Order;
   //Finally we need to renumber the priorities if the ring doesn't have
   //symmetry
   Pri_t::iterator Pi=Priors.begin();
   Pri_t::reverse_iterator rPi=Priors.rbegin();
   ++Pi;//Need to move to the two position
   size_t NComps=(Priors.size()-Priors.size()%2)/2,counter=0;
   for(;counter<NComps;++counter,++Pi,++rPi)
      if(*Pi!=*rPi)break;
   Symm=true;
   if(counter<NComps){
      Symm=false;
      Pri_t Temp;
      for(counter=0,Pi=Priors.begin();Pi!=Priors.end();++Pi,++counter)
         Temp.push_back(counter);
      Priors=Temp;
   }
   return Priors;
}




}}//end namespaces




#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_ */
