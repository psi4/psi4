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
#include <algorithm>
#include "Node.h"
#include "NodeUtils.h"
#include "Exception2.h"

namespace psi {
namespace LibMolecule {
typedef boost::shared_ptr<Node> SharedNode;
typedef std::vector<SharedNode> vSN_t;
typedef vSN_t::iterator vSN_Itr;

//Function to compare two nodes (just compares the addresses)
static bool Compare(const SharedNode LHS, const SharedNode RHS) {
   return LHS.get()<RHS.get();
}

SharedNode Node::GetConnSubNode(SharedNode Other) {
   ConnItr It2=ConnBegin(),It2End=ConnEnd();
   for (;It2!=It2End; ++It2)
      if(It2->second.get()==Other.get()) return It2->first;
   return Other;
}

unsigned Node::Z() const {
   unsigned Max=0;
   Node::const_iterator It=SubBegin(),ItEnd=SubEnd();
   for (; It!=ItEnd; ++It)Max=std::max((*It)->Z(),Max);
   return Max;
}

unsigned Node::size() const {
   return MyMember_.size();
}

int Node::operator[](const unsigned i) const {
   return MyMember_[i];
}

void Node::AddConn(SharedNode NewNode, SharedNode OldNode) {
   if (OldNode.get()==this)
      ConnNodes_.push_back(Pair_t(OldNode,NewNode));
   else {
      ConnItr It=ConnBegin(),ItEnd=ConnEnd();
      for (;It!=ItEnd; ++It)
         if (It->second.get()==OldNode.get()) {
            It->second=NewNode;
            break;
         }
   }
}

void Node::AddTypes(const ParamT& Parent,
     const PsiMap<size_t, size_t>& Order,
     size_t& Priority,
     bool Symm,
     SharedNode_t Left) {
   size_t TempPrior=Priority;
   Node::iterator It=SubBegin(),ItEnd=SubEnd();
   if (It==ItEnd){//I am a primitive node...
      MyTypes_.push_back(
            ParamT(MyTypes_.back(),
                  Parent.Base()[0],
                  Parent.Base()[1],
                  Order,Priority));
   }
   size_t TempPriority=Priority,MaxPriority=0,Offset=0;
   for(size_t counter=0;It!=ItEnd;++It,++counter){
      Offset=Prior_[counter]+(IsSymm_&&!Symm&&!IsPrim()?counter:0);
      TempPriority+=Offset;
      (*It)->AddTypes(Parent,Order,TempPriority,Symm,Left);
      MaxPriority=(TempPriority>MaxPriority?TempPriority:MaxPriority);
      TempPriority-=Offset;
   }
   TempPrior+=Offset;
   Priority=(IsSymm_&&!IsPrim()?TempPrior:MaxPriority);
}

std::vector<size_t> Node::DeterminePrior(
      std::deque<SharedNode>& FN,
      bool& Symm) const {
   size_t size=FN.size();
   std::vector<size_t> Priors(size, 0);
   typedef std::deque<SharedNode> FN_t;
   FN_t::const_iterator It=FN.begin(),ItEnd=FN.end();
   FN_t::const_reverse_iterator rIt=FN.rbegin(),rItEnd=FN.rend();
   size_t NComp=(size-size%2)/2;
   size_t CurrPrior=0;
   bool Reorder=false,LLessR=false;
   Symm=true;
   //Compare them pairwise
   for (size_t i=0; i<NComp; ++i, ++It, ++rIt) {
      size_t Lidx=i,Ridx=size-1-i;
      //As soon as orderchecked is true we have broken the symmetry
      if ((*It)->Type()==(*rIt)->Type())
         Priors[Lidx]=Priors[Ridx]=CurrPrior++;
      else {
         Symm=false;
         const Node* Result=Priority(It->get(),rIt->get());
         if(Result==NULL)PSIERROR("Couldn't resolve priority");
         LLessR=(Result==It->get());
         if(!LLessR) Reorder=true;
         for(size_t j=0;j<Priors.size();j++)Priors[j]=j;
         break;
      }
   }
   if(size%2==1&&Symm) Priors[NComp]=CurrPrior;
   if(Reorder){
      //std::vector<size_t> TempPriors(Priors.rbegin(),Priors.rend());
      //Priors=TempPriors;
      FN_t temp(FN.rbegin(),FN.rend());
      FN=temp;
   }
   return Priors;
}

template<typename T>
static void Flip(T& In){
   T Temp(In.rbegin(),In.rend());
   In=Temp;
}

void Node::FillNode(std::deque<SharedNode>& FoundNodes) {
   typedef std::deque<SharedNode> FN_t;
   typedef FN_t::const_iterator FN_tItr;
   IsSymm_=true;
   Prior_=DeterminePrior(FoundNodes,IsSymm_);
   ConnNodes_=DetermineConns(FoundNodes.begin(),FoundNodes.end());
   for(FN_tItr It=FoundNodes.begin(); It!=FoundNodes.end(); ++It)AddSubNode(*It);
   ParamT orig=MyTypes_.back();
   MyTypes_.pop_back();
   Node::ConnItr It1=ConnBegin(),It1End=ConnEnd();
   PsiMap<size_t, size_t> Orders;
   for (; It1!=It1End; ++It1) {
      Node::iterator It2=SubBegin(),It2End=SubEnd();
      for (size_t Counter=0; It2!=It2End; ++It2,++Counter)
         if(It1->first.get()==It2->get()){
               ActiveSubNodes_.push_back(Counter);
               Orders[Counter]++;
         }
   }
   std::string Prefix="",LongPrefix="";
   if((*SubNodes_.begin())->Type().Base()[0]=="CT"){
      Prefix="CT";
      LongPrefix="C-Terminal ";
   }
   else if((*SubNodes_.begin())->Type().Base()[0]=="NT"){
      Prefix="NT";
      LongPrefix="N-Terminal ";
   }
   MyTypes_.push_back(ParamT(
         Prefix+orig.Base()[0],
         LongPrefix+orig.Base()[1], Orders));
   std::vector<size_t>::iterator Pi,PiEnd=Prior_.end();
   if(Orders.size()>1&&IsSymm_){
      //If we have more than one connection point and the priorities are
      //symmetric we have an orientation ambiguity
      size_t NComps=(Prior_.size()-Prior_.size()%2)/2;
      bool reorder=false;
      ConnItr ConnI=ConnNodes_.begin();
      std::vector<Pair_t>::reverse_iterator rConnI=ConnNodes_.rbegin();
      NComps=(ConnNodes_.size()-ConnNodes_.size()%2)/2;
      for(size_t counter=0;counter<NComps&&!reorder;
            ++counter,++ConnI,++rConnI){
         const Node* Result=Priority(ConnI->second.get(),
                                    rConnI->second.get());
         if(Result==NULL)continue;
         else if(Result==ConnI->second.get())break;
         else reorder=true;
      }
      if(reorder){
         //For good measure flip (almost) everything
         Flip(Prior_);
         //We don't flip this or else a 0,1 double bond shows up as 1,0
         //which due to symmetry, is just a bad numbering...
         //Flip(ActiveSubNodes_);
         Flip(ConnNodes_);
         Flip(SubNodes_);
         Flip(PrimNodes_);
      }
   }
   Pi=Prior_.begin();PiEnd=Prior_.end();
   size_t TempPrior=*Pi;
   Node::iterator It2=SubBegin(),It2End=SubEnd();
   SharedNode_t LastNode;
   for (; It2!=It2End; ++It2,++Pi) {
      TempPrior=(IsSymm_?*Pi:TempPrior);
      (*It2)->AddTypes(MyTypes_.back(), Orders,TempPrior,IsSymm_,LastNode);
      TempPrior+=1;
      LastNode=*It2;
   }
}

Node::Node(const std::string& BaseAbbrv, const std::string& BaseName) :
      MyTypes_(1, ParamT(BaseAbbrv, BaseName)),Prior_(1,0),IsSymm_(true) {
}

void Node::AddSubNode(boost::shared_ptr<Node> NewNode) {
   SubNodes_.push_back(NewNode);
   if (NewNode->PrimNodes_.size()>0)
      PrimNodes_.insert(PrimNodes_.end(),NewNode->PrimBegin(), NewNode->PrimEnd());
   else
      PrimNodes_.push_back(NewNode);

   //Add the members of our new subnode to the array
   for (size_t i=0; i<SubNodes_.back()->size(); i++)
      MyMember_.push_back((*SubNodes_.back())[i]);
}


std::string Node::PrintOut(const std::string spaces, const size_t Level) const {
   std::stringstream Message;
   std::string MyStrType=MyTypes_.back().PrintOut();
   std::string NewSpaces="    ";
   //Print my atoms
   if (spaces==""||SubNodes_.size()==0) {
      Message<<spaces;
      /* This is taken care of in ParameterType.cc
        if (ActiveSubNodes_.size()>1) {
         std::vector<size_t>::const_iterator
         It=ActiveSubNodes_.begin(),ItEnd=ActiveSubNodes_.end();
         --ItEnd;
         for (; It!=ItEnd; ++It)
            Message<<(*It)<<",";
         Message<<(*It)<<"-";
      }*/
      Message<<MyStrType<<" Atom(s):";
      for (unsigned MemberI=0; MemberI<size(); MemberI++)
         Message<<(*this)[MemberI]<<" ";
      Message<<std::endl;
   }
   if (Level==0) return Message.str();
   //Print my subnodes
   vSN_t::const_iterator It=SubNodes_.begin(),ItEnd=SubNodes_.end();
   for (; It!=ItEnd; ++It)
      Message<<(*It)->PrintOut(NewSpaces);
   //Print my type hierarchy
   if (SubNodes_.size()==0) {
      int NTypes=MyTypes_.size();
      if (NTypes>1) Message<<NewSpaces<<MyTypes_[0].PrintOut()<<" --> ";
      for (int TypeI=1; TypeI<NTypes-1; TypeI++)
         Message<<MyTypes_[TypeI].PrintOut()<<" --> ";
      if (NTypes>1) Message<<MyTypes_[NTypes-1].PrintOut()<<std::endl;
   }
   //If I'm the top-level node print my connections
   if (spaces=="") {
      Message<<"My Connections: "<<std::endl;

      const_ConnItr It1=ConnBegin(),It1End=ConnEnd();
      for (; It1!=It1End; ++It1) {
         Message<<NewSpaces<<It1->second->Type().PrintOut()<<" Atom(s):";
         for (unsigned MemberI=0; MemberI<It1->second->size(); MemberI++)
            Message<<(*(It1->second))[MemberI]<<" ";
         Message<<std::endl;
      }
   }
   return Message.str();
}


void Node::UpdateConns(SharedNode NewMe) {
   ConnItr It=ConnBegin(),ItEnd=ConnEnd();
   Node::iterator It2,It2End=SubEnd();
   for (; It!=ItEnd; ++It)
      for (It2=SubNodes_.begin(); It2!=It2End; ++It2)
         It->second->AddConn(NewMe, (*It2));
}

}
} //End namespaces
