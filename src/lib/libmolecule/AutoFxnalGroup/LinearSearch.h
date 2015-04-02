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
      virtual FxnGrpType LHSType(const unsigned)const{return FxnGrpType::NONE;}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N+1-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual FxnGrpType RHSType(const unsigned)const{return FxnGrpType::NONE;}
      boost::shared_ptr<Node> Reorder(boost::shared_ptr<Node> NodeIn);
   public:
      ///Same arguments as Node's constructor
      template<typename...MyTypes>
      LinearSearch<CurrentType,LinearTypes...>
         (MyTypes...Types):Base_t(Types...){
      }
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
      template<typename...MyTypes>
      LinearSearch<LastType>
         (MyTypes...Types):Node(Types...){
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
      virtual FxnGrpType LHSType(const unsigned)const{return FxnGrpType::NONE;}
      /** \brief Function for specifying the subtype of the connection
       *
       *   Assume our queue has N elements and we attempting to add the
       *   N+1-th element to it.  This means we have found N-1 connections
       *   and are attempting to add the N-th connection.  Given N, this
       *   function should be overridden to return the type the N+1-th
       *   element's connecting subnode must have.  By default this is
       *   NONE, which means the check is ignored.
       */
      virtual FxnGrpType RHSType(const unsigned)const{return FxnGrpType::NONE;}
      boost::shared_ptr<Node> Reorder(boost::shared_ptr<Node> NodeIn);
   public:
      ///Trivial check of Node vs. LastType
      virtual bool FindMe(boost::shared_ptr<Node> Node);
};

/*******************Implementations********************************/
template<typename CurrentType,typename...LinearTypes>
boost::shared_ptr<Node> LinearSearch<CurrentType,LinearTypes...>::
Reorder(boost::shared_ptr<Node> NodeIn){
   boost::shared_ptr<CurrentType> Temp(new CurrentType());
   if(Temp.FindMe(NodeIn))
      return Temp;
   else return boost::shared_ptr<CurrentType>();
}

template<typename CurrentType>
boost::shared_ptr<Node> LinearSearch<CurrentType>::
Reorder(boost::shared_ptr<Node> NodeIn){
   boost::shared_ptr<CurrentType> Temp(new CurrentType());
   if(Temp.FindMe(NodeIn))
      return Temp;
   else return boost::shared_ptr<CurrentType>();
}

template<typename CurrentType,typename...LinearTypes>
bool LinearSearch<CurrentType,LinearTypes...>::FindMe(boost::shared_ptr<Node> Nodes){
   typedef boost::shared_ptr<Node> SharedNode;
   CurrentType Temp;
   if(Temp.Type()!=Nodes->Type())return false;
   std::deque<SharedNode> NodeQueue({Nodes});
   if(Base_t::FindMe(NodeQueue)){
      return true;
   }
   return false;
}

static bool IsGood(
      boost::shared_ptr<Node> NodeI,
      boost::shared_ptr<Node> NodeJ,const Node& Temp,
      const std::deque<boost::shared_ptr<Node> >& FNs,
      const FxnGrpType& LHS,
      const FxnGrpType& RHS){
   typedef boost::shared_ptr<Node> SharedNode;
   if(NodeJ->Type()!=Temp.Type())return false;
   std::deque<SharedNode>::const_iterator GI=FNs.begin(),GIEnd=FNs.end();
   for(;GI!=GIEnd;++GI)
      if(GI->get()==NodeJ.get())return false;
   if(LHS!=FxnGrpType::NONE){
      SharedNode Conn=NodeI->GetConnSubNode(NodeJ);
      if(Conn.get()==NodeJ.get()||Conn->Type()!=LHS){
         return false;
      }
   }
   if(RHS!=FxnGrpType::NONE){
      SharedNode Conn=NodeJ->GetConnSubNode(NodeI);
      if(Conn.get()==NodeI.get()||Conn->Type()!=RHS){
         return false;
      }
   }
   return true;
}

template<typename CurrentType, typename...LinearTypes>
bool LinearSearch<CurrentType,LinearTypes...>::
   FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes){
   typedef boost::shared_ptr<Node> SharedNode;
   typedef std::vector<SharedNode> Conn_t;
   SharedNode NodeI=(IsLinear()?FoundNodes.back():FoundNodes.front());
   CurrentType Temp;
   Conn_t::iterator It=NodeI->ConnNodes_.begin(),ItEnd=NodeI->ConnNodes_.end();
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
   if(Temp.Type()==FoundNodes->Type())this->AddSubNode(FoundNodes);
   else return false;
   return true;
}
template<typename LastType>
bool LinearSearch<LastType>::
   FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes){
   typedef boost::shared_ptr<Node> SharedNode;
   typedef std::vector<SharedNode> Conn_t;
   SharedNode NodeI=(IsLinear()?FoundNodes.back():FoundNodes.front());
   LastType Temp;
   Conn_t::iterator It=NodeI->ConnNodes_.begin(),ItEnd=NodeI->ConnNodes_.end();
   for(;It!=ItEnd;++It){
      if(!IsGood(NodeI,*It,Temp,FoundNodes,
            LHSType(FoundNodes.size()),
            RHSType(FoundNodes.size())))continue;
      FoundNodes.push_back(*It);
      while(!FoundNodes.empty()){
         this->AddSubNode(FoundNodes.front());
         FoundNodes.pop_front();
      }
      return true;
   }
   return false;
}

}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_LINEARSEARCH_H_ */
