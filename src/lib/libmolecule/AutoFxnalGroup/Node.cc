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
#include "exception.h"


namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<Node> SharedNode;
typedef std::vector<SharedNode> vSN_t;
typedef vSN_t::iterator vSN_Itr;

//Function to compare two nodes (just compares the addresses)
static bool Compare(const SharedNode LHS,const SharedNode RHS){
   return LHS.get()<RHS.get();
}

SharedNode Node::GetConnSubNode(SharedNode Other){
   vSN_Itr It=SubNodes_.begin(),ItEnd=SubNodes_.end();
   for(;It!=ItEnd;++It){
      vSN_Itr It2=(*It)->ConnNodes_.begin(),It2End=(*It)->ConnNodes_.end();
      for(;It2!=It2End;++It2)
         if(It2->get()==Other.get())return *It;
   }
   return Other;
}

unsigned Node::Z()const{
   unsigned max=0;
   vSN_t::const_iterator It=SubNodes_.begin(),ItEnd=SubNodes_.end();
   for(;It!=ItEnd;++It)max=((*It)->Z()>max?(*It)->Z():max);
   return max;
}

unsigned Node::size()const{return MyMember_.size();}

int Node::operator[](const unsigned i)const{return MyMember_[i];}

void Node::AddConn(SharedNode NewNode, SharedNode OldNode){
   if(!OldNode)ConnNodes_.push_back(NewNode);
   else{
      vSN_Itr It=ConnNodes_.begin(),ItEnd=ConnNodes_.end();
      for(;It!=ItEnd;++It)
         if(It->get()==OldNode.get()){
            (*It)=NewNode;
            break;
         }
   }
   std::sort(ConnNodes_.begin(),ConnNodes_.end(),Compare);
}

void Node::AddTypes(const ParamT& Parent,
                    const size_t Order,
                    const size_t Priority){
   if(SubNodes_.size()==0){
      MyTypes_.push_back(ParamT(MyTypes_.front(),
                               Parent.Base()[0],
                               Parent.Base()[1],
                               Order,Priority));
   }
   else{
      vSN_Itr SubNodeI=SubNodes_.begin(),SubNodeIEnd=SubNodes_.end();
      for(;SubNodeI!=SubNodeIEnd;++SubNodeI)
         (*SubNodeI)->AddTypes(Parent,Order,Priority);
   }
}

static void Union(vSN_t& LHS, vSN_t&RHS, vSN_t& Result){
   vSN_Itr It=std::set_union(LHS.begin(),LHS.end(),
         RHS.begin(),RHS.end(),Result.begin(),Compare);
   Result.resize(It-Result.begin());
}

void Node::AddSubNode(boost::shared_ptr<Node> NewNode,
                      const size_t Order,
                      const size_t Prior){
   SubNodes_.push_back(NewNode);
   if(NewNode->PrimNodes_.size()>0){
      vSN_t TempPrims(PrimNodes_.size()+NewNode->PrimNodes_.size());
      Union(PrimNodes_,NewNode->PrimNodes_,TempPrims);
      PrimNodes_=TempPrims;
   }
   else PrimNodes_.push_back(NewNode);
   SubNodes_.back()->AddTypes(MyTypes_.back(),Order,Prior);
   //Add the members of our new subnode to the array
   for(size_t i=0;i<SubNodes_.back()->size();i++)
      MyMember_.push_back((*SubNodes_.back())[i]);
   //Get a list of the nodes we are now connected to
   vSN_t& OtherNodes=NewNode->ConnNodes_;
   vSN_t Result(ConnNodes_.size()+OtherNodes.size());
   Union(ConnNodes_,OtherNodes,Result);
   //Remove the connections that are part of this

   //Don't want to sort actual SubNodes array as it destroys the order
   vSN_t Temp=SubNodes_,Result2(Result.size());
   std::sort(Temp.begin(),Temp.end(),Compare);
   vSN_Itr It=std::set_difference(Result.begin(),Result.end(),
         Temp.begin(),Temp.end(),Result2.begin(),Compare);
   Result2.resize(It-Result2.begin());
   ConnNodes_=Result2;
}

std::string Node::PrintOut(const std::string spaces)const{
   std::stringstream Message;
   std::string MyStrType=MyTypes_.back().PrintOut();
   std::string NewSpaces="    ";
   //Print my atoms
   if(spaces==""||SubNodes_.size()==0){
      Message<<spaces;
      if(ActiveSubNodes_.size()>1){
         std::set<size_t>::const_iterator It=ActiveSubNodes_.begin(),
               ItEnd=ActiveSubNodes_.end();
         --ItEnd;
         for(;It!=ItEnd;++It)
            Message<<(*It)<<",";
         Message<<(*It)<<"-";
      }
      Message<<MyStrType<<" Atom(s):";
      for(unsigned MemberI=0;MemberI<size();MemberI++)
         Message<<(*this)[MemberI]<<" ";
      Message<<std::endl;
   }
   //Print my subnodes
   vSN_t::const_iterator It=SubNodes_.begin(),ItEnd=SubNodes_.end();
   for(;It!=ItEnd;++It) Message<<(*It)->PrintOut(NewSpaces);
   //Print my type hierarchy
   if(SubNodes_.size()==0){
      int NTypes=MyTypes_.size();
      if(NTypes>1)
         Message<<NewSpaces<<MyTypes_[0].PrintOut()<<" --> ";
      for(int TypeI=1;TypeI<NTypes-1;TypeI++)
         Message<<MyTypes_[TypeI].PrintOut()<<" --> ";
      if(NTypes>1)
         Message<<MyTypes_[NTypes-1].PrintOut()
                 <<std::endl;
   }
   //If I'm the top-level node print my connections
   if(spaces==""){
      Message<<"My Connections: "<<std::endl;

      It=ConnNodes_.begin();ItEnd=ConnNodes_.end();
      for(;It!=ItEnd;++It){
         Message<<NewSpaces<<(*It)->Type().PrintOut()<<" Atom(s):";
         for(unsigned MemberI=0;MemberI<(*It)->size();MemberI++)
            Message<<(*(*It))[MemberI]<<" ";
         Message<<std::endl;
      }
   }
   return Message.str();
}

void Node::UpdateConns(SharedNode NewMe){
   vSN_Itr It=ConnNodes_.begin(),ItEnd=ConnNodes_.end(),It2,
         It2End=SubNodes_.end();
   for(;It!=ItEnd;++It){
      for(It2=SubNodes_.begin();It2!=It2End;++It2)
         (*It)->AddConn(NewMe,(*It2));
   }
}

}}//End namespaces
