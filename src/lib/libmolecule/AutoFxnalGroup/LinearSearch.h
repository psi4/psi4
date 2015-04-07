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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_LINEARSEARCH_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_LINEARSEARCH_H_
#include <deque>
#include "Node.h"
namespace psi{
namespace LibMolecule{

template<typename CurrentType, typename... LinearTypes>
class LinearSearch : public LinearSearch<LinearTypes...>{
   protected:
      ///Convenient typedef of the base class
      typedef LinearSearch<LinearTypes...> Base_t;
      ///The internal call that continues the recursion
      bool FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes);
      ///A function used to flag if this is a radial or linear search
      virtual bool IsLinear()const{return true;}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual ParamT LHSType(const unsigned)const{
         return ParamT("NONE","NONE");}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N+1-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual ParamT RHSType(const unsigned)const{
         return ParamT("NONE","NONE");}
   public:
      ///Same arguments as Node's constructor
      LinearSearch<CurrentType,LinearTypes...>
         (const std::string& BaseAbbv,const std::string& BaseName):
         Base_t(BaseAbbv,BaseName){}
      ///STFU compiler...
      virtual ~LinearSearch<CurrentType,LinearTypes...>(){}
      ///Given a starting node, tries to find the desired group
      virtual bool FindMe(boost::shared_ptr<Node> Nodes);
};

template<typename LastType>
class LinearSearch<LastType> : public Node{
   protected:
      ///The internal call that continues the recursion
      bool FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes);

      virtual bool IsLinear()const{return true;}
      ///Same arguments as Node's constructor
      LinearSearch<LastType>
         (const std::string& Abbv,const std::string& Name):Node(Abbv,Name){
      }
      ///STFU compiler...
      virtual ~LinearSearch<LastType>(){}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual ParamT LHSType(const unsigned)const{return ParamT("NONE","NONE");}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N+1-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual ParamT RHSType(const unsigned)const{return ParamT("NONE","NONE");}
   public:
      ///Trivial check of Node vs. LastType
      virtual bool FindMe(boost::shared_ptr<Node> Node);
};

/*******************Implementations********************************/
static bool IsGood(
      boost::shared_ptr<Node> NodeI,
      boost::shared_ptr<Node> NodeJ,const Node& Temp,
      const std::deque<boost::shared_ptr<Node> >& FNs,
      const ParamT& LHS,
      const ParamT& RHS){
   typedef boost::shared_ptr<Node> SharedNode;
   size_t i=0;
   for(;i<Temp.NTypes();i++)
      if(NodeJ->Type()==Temp.Type(i))break;
   if(i==Temp.NTypes())return false;
   std::deque<SharedNode>::const_iterator GI=FNs.begin(),GIEnd=FNs.end();
   for(;GI!=GIEnd;++GI)
      if(GI->get()==NodeJ.get())return false;
   if(LHS!=ParamT("NONE","NONE")){
      SharedNode Conn=NodeI->GetConnSubNode(NodeJ);
      if(Conn.get()==NodeJ.get()||Conn->Type()!=LHS){
         return false;
      }
   }
   if(RHS!=ParamT("NONE","NONE")){
      SharedNode Conn=NodeJ->GetConnSubNode(NodeI);
      if(Conn.get()==NodeI.get()||Conn->Type()!=RHS){
         return false;
      }
   }
   return true;
}

template<typename CurrentType,typename...LinearTypes>
bool LinearSearch<CurrentType,LinearTypes...>::FindMe(boost::shared_ptr<Node> Nodes){
   typedef boost::shared_ptr<Node> SharedNode;
   CurrentType Temp;
   std::deque<SharedNode> EmptyQueue;
   if(!IsGood(Nodes,Nodes,Temp,EmptyQueue,
         LHSType(EmptyQueue.size()),
         RHSType(EmptyQueue.size())))return false;
   std::deque<SharedNode> NodeQueue({Nodes});
   if(Base_t::FindMe(NodeQueue)){
      return true;
   }
   return false;
}



template<typename CurrentType, typename...LinearTypes>
bool LinearSearch<CurrentType,LinearTypes...>::
   FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes){
   typedef boost::shared_ptr<Node> SharedNode;
   typedef std::vector<SharedNode> Conn_t;
   SharedNode NodeI=(IsLinear()?FoundNodes.back():FoundNodes.front());
   CurrentType Temp;
   Node::iterator It=NodeI->ConnBegin(),ItEnd=NodeI->ConnEnd();
   for(;It!=ItEnd;++It){
      if(!IsGood(NodeI,*It,Temp,FoundNodes,
                 LHSType(FoundNodes.size()),
                 RHSType(FoundNodes.size())))continue;
      FoundNodes.push_back(*It);
      if(!Base_t::FindMe(FoundNodes))FoundNodes.pop_back();
      else return true;
   }
   return false;
}
template<typename LastType>
bool LinearSearch<LastType>::
   FindMe(boost::shared_ptr<Node> FoundNodes){
   LastType Temp;
   if(Temp.Type()==FoundNodes->Type()){
      size_t Order=FoundNodes->NEdges();
      this->AddSubNode(FoundNodes,Order,0);
      ParamT Orig=MyTypes_[0];
      this->MyTypes_[0]=ParamT(Orig.Base()[0],Orig.Base()[1],Order,0);
   }
   else return false;
   return true;
}

static void PriorOrder(
      const std::deque<boost::shared_ptr<Node> >& FoundNodes,
      int& Prior, size_t& Order,
      std::set<size_t>& ActiveSubNodes){
   std::deque<boost::shared_ptr<Node> >::const_iterator
         It1=FoundNodes.begin(),It2,
         It1End=FoundNodes.end();
   std::set<size_t> AllNodes;
   for(size_t counter=0;It1!=It1End;++It1,++counter){
      Node::iterator ItEnd=(*It1)->ConnEnd();
      Node::iterator It=(*It1)->ConnBegin();
      for(;It!=ItEnd;++It){
        It2=FoundNodes.begin();
        for(;It2!=It1End;++It2)
           if(It2->get()==It->get())break;
        if(It2!=It1End)continue;
        AllNodes.insert(counter);
        Order++;
     }
   }
   ActiveSubNodes=AllNodes;
}

template<typename LastType>
bool LinearSearch<LastType>::
   FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes){
   typedef boost::shared_ptr<Node> SharedNode;
   SharedNode NodeI=(IsLinear()?FoundNodes.back():FoundNodes.front());
   LastType Temp;
   Node::iterator It=NodeI->ConnBegin(),ItEnd=NodeI->ConnEnd();
   for(;It!=ItEnd;++It){
      if(!IsGood(NodeI,*It,Temp,FoundNodes,
            LHSType(FoundNodes.size()),
            RHSType(FoundNodes.size())))continue;
      FoundNodes.push_back(*It);
      int Prior=-1;
      size_t Order=0;
      PriorOrder(FoundNodes,Prior,Order,this->ActiveSubNodes_);
      ParamT Front=this->Type();
      this->MyTypes_[MyTypes_.size()-1]=
            ParamT(Front,Front.Base()[0],Front.Base()[1],Order,-1);
      while(!FoundNodes.empty()){
         this->AddSubNode(FoundNodes.front(),Order,Prior);
         FoundNodes.pop_front();
      }
      return true;
   }
   return false;
}

}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_LINEARSEARCH_H_ */
