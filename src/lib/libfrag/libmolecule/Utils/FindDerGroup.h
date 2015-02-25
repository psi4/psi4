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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FINDDERGROUP_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FINDDERGROUP_H_
#include<vector>
#include "OrganicGeom.h"
namespace psi{
namespace LibMolecule{

///Class to hold code factorization of the forthcoming classes
template<typename GroupName,typename T>
class FindDerGroupHelp{
   protected:
      ///Returns true if we found a group with type T bonded to atom "index"
      static bool CheckGroup(int& index,
            std::vector<int>& Groups,const Connections& Conns,
            ConnGroups& FoundGroups,bool Prim=false){
         bool IsGood=false;
         for(int i=0;i<Conns[index].size();i++){
            int Conn=Conns[index][i];
            //Various checks for:
            //Already found
            bool found=false;
            for(int j=0;j<Groups.size()&&!found;j++){
            	if(Groups[j]==Conn)found=true;
            }
            if(found)continue;
            //Already part of a functional group
            if(FoundGroups.count(Conn)==0)continue;
            //Is literally the same group as our index
            if(&(*FoundGroups[Conn])==&(*FoundGroups[index]))continue;
            //Is the group we want
            if(FoundGroups[Conn]->Type()==T().Type()){
               Groups.push_back(Conn);
               IsGood=true;
               if(!Prim)index=Conn;
               break;
            }
         }
         return IsGood;
      }
      ///Takes Groups, makes a "GroupName", and updates "FoundGroups"
      static void UpdateGroups(std::vector<int>& Groups,ConnGroups& FoundGroups){
         std::vector<boost::shared_ptr<const FxnalGroup> >
            ActualGroups(Groups.size());
         std::vector<int>::const_iterator GroupI=Groups.begin(),
                                     GroupEnd=Groups.end();
         for(int i=0;GroupI!=GroupEnd;++GroupI,i++){
            ActualGroups[i]=FoundGroups[(*GroupI)];
            FoundGroups.erase(*GroupI);
         }
         boost::shared_ptr<GroupName> temp(new GroupName());
         temp->SetGroups(ActualGroups);
         for(int i=0;i<temp->NAttachPoint();i++)
            FoundGroups[temp->AttachPoint(i)]=temp;
         if(temp->NAttachPoint()==0){
            FoundGroups[(*Groups.begin())]=temp;
         }
      }
};

/*

template<typename T>
class FindDerGroupHelp<Carbonyl,T>{
   protected:
      typedef Carbonyl GroupName;
      ///Returns true if we found the group
      static bool CheckGroup(int& index,
            std::set<int>& Groups,const Connections& Conns,
            ConnGroups& FoundGroups,bool Prim=false){
         bool IsGood=false;
         for(int i=0;i<Conns[index].size();i++){
            int Conn=Conns[index][i];
            std::cout<<"Checking connection: "<<Conn<<std::endl;
            //Various checks for:
            //Already found
            if(Groups.count(Conn)==1)continue;
            std::cout<<"Not Found already"<<std::endl;
            //Already part of a functional group
            if(FoundGroups.count(Conn)==0)continue;
            std::cout<<"And it's an attachment point still"<<std::endl;
            //Is literally the same group as our index
            if(&(*FoundGroups[Conn])==&(*FoundGroups[index]))continue;
            std::cout<<"Not literally the same group"<<std::endl;
            //Is the group we want
            if(FoundGroups[Conn]->Type()==T().Type()){
               std::cout<<"It's the type we want"<<std::endl;
               Groups.insert(Conn);
               IsGood=true;
               if(!Prim)index=Conn;
               break;
            }
         }
         return IsGood;
      }

      static void UpdateGroups(std::set<int>& Groups,ConnGroups& FoundGroups){
         std::vector<boost::shared_ptr<const FxnalGroup> >
            ActualGroups(Groups.size());
         std::set<int>::const_iterator GroupI=Groups.begin(),
                                     GroupEnd=Groups.end();
         for(int i=0;GroupI!=GroupEnd;++GroupI,i++){
            ActualGroups[i]=FoundGroups[(*GroupI)];
            FoundGroups.erase(*GroupI);
         }
         boost::shared_ptr<GroupName> temp(new GroupName());
         temp->SetGroups(ActualGroups);
         for(int i=0;i<temp->NAttachPoint();i++)
            FoundGroups[temp->AttachPoint(i)]=temp;
         if(temp->NAttachPoint()==0){
            FoundGroups[(*Groups.begin())]=temp;
         }
      }
};
*/

/** \brief A functor that can assign derived groups
 *
 *  We define two types of derived groups: primitive groups, and what
 *  we will just call derived groups
 *
 *  For our purposes a primitive group is an element and a variable number
 *  of hydrogens.  More generally a primitive group is some group X, and
 *  and a variable number of other groups \f$Y_i\f$, such that all \f$Y_i\f$
 *  are attached to X.  Hence for each type in the parameter pack, Y, we see
 *  if Y is attached to X and if so add it to a growing array of found groups.
 *
 *  For example a Methyl looks like Carbon, Hydrogen,
 *  Hydrogen, Hydrogen, Hydrogen.  The heteroatom should always be first
 *  the first type passed to this functor (after the name of type of the
 *  resulting group).
 *
 *  Derived groups are
 *
 *
 *
 */
template<typename GroupName,
    typename ActiveGroup, typename...Args>
class FindDerGroup:
      public FindDerGroup<GroupName,Args...>,
      private virtual FindDerGroupHelp<GroupName,ActiveGroup>{
   private:
      typedef FindDerGroup<GroupName,Args...> Base1_t;
      typedef FindDerGroupHelp<GroupName,ActiveGroup> Base2_t;
   protected:
      static bool FindGroup(int index,std::vector<int>& Groups,
                      const Connections& Conns,ConnGroups& FoundGroups,
                      bool Prim=false){
         if(!Base2_t::CheckGroup(index,Groups,Conns,FoundGroups,Prim))
            return false;
         return Base1_t::FindGroup(index,Groups,Conns,FoundGroups,Prim);
      }
   public:
      static bool FindGroup(boost::shared_ptr<const FxnalGroup> Group,
            const Connections& Conns,ConnGroups& FoundGroups,bool Prim=false){
         int NAttach=Group->NAttachPoint();
         ActiveGroup Temp;
         if(Group->Type()!=Temp.Type())return false;
         bool found=false;
         for(int i=0;i<NAttach;i++){
            int Index=Group->AttachPoint(i);
            std::vector<int> Groups;
            if(Prim&&Conns[Index].size()<sizeof...(Args))continue;
            Groups.push_back(Index);
            found=Base1_t::FindGroup(Index,Groups,Conns,FoundGroups,Prim);
            if(found)break;
         }
         return found;
      }
};

///Recursion end point
template<typename GroupName, typename ActiveGroup>
class FindDerGroup<GroupName,ActiveGroup>:
      private virtual FindDerGroupHelp<GroupName,ActiveGroup>{
   private:
      typedef FindDerGroupHelp<GroupName,ActiveGroup> Base_t;
   protected:
      static bool FindGroup(int index,std::vector<int>& Groups,
            const Connections& Conns,ConnGroups& FoundGroups,bool Prim=false){
         if(!Base_t::CheckGroup(index,Groups,Conns,FoundGroups,Prim))return false;
         Base_t::UpdateGroups(Groups,FoundGroups);
         return true;
      }
};

template<int Valency,typename GroupName,typename...Args>
class FindPrimGroup:public FindDerGroup<GroupName,Args...>{
   private:
      typedef FindDerGroup<GroupName,Args...> Base;
   public:
      static  bool FindGroup(boost::shared_ptr<const FxnalGroup> Group,
            const Connections& Conns,ConnGroups& FoundGroups){
         if(Conns[Group->AttachPoint()].size()!=Valency)return false;
         return Base::FindGroup(Group,Conns,FoundGroups,true);
      }
};

template<int Valency,typename GroupName,typename T>
class FindPrimGroup<Valency,GroupName,T>:
   public FindDerGroup<GroupName,T>,
   private virtual FindDerGroupHelp<GroupName,T>{
   private:
      typedef FindDerGroup<GroupName,T> Base;
      typedef FindDerGroupHelp<GroupName,T> Base2;
   public:
      static  bool FindGroup(boost::shared_ptr<const FxnalGroup> Group,
            const Connections& Conns,ConnGroups& FoundGroups){
         if(Conns[Group->AttachPoint()].size()!=Valency)return false;
         if(Group->Type()!=T().Type())return false;
         std::vector<int> temp;
         temp.push_back(Group->AttachPoint());
         Base2::UpdateGroups(temp,FoundGroups);
         return true;
      }
};


}}//End namespace



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FINDDERGROUP_H_ */
