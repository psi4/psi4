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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDGROUPFUNCTOR_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDGROUPFUNCTOR_H_
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "FxnalGroup.h"
#include "PsiMap.h"
namespace psi{
namespace LibMolecule{
typedef PsiMap<int,boost::shared_ptr<const FxnalGroup> > PsiMap_t;
typedef std::vector<std::vector<int> > Conn_t;
typedef boost::shared_ptr<const Molecule> SharedMol;
template<typename T>
void UpdateGroups(PsiMap_t& List, T Object){
   List[Object[0]]=boost::shared_ptr<T>(new T(Object));
}


class NullGroup{
   public:
      NullGroup(int*){}
      FxnGroup_t Type()const{return NO_GROUP;}
};
template<>
void UpdateGroups<NullGroup>(PsiMap_t&,NullGroup){
   LibMoleculeBase BSClass;
      BSClass.Error("The atom you choose isn't allowed to make this many "
                    "connections");
}

template<int N,typename Base,typename T1,typename T1DB=NullGroup,typename T2T=NullGroup,
         typename T2=NullGroup,typename T2DB=NullGroup,typename T2TB=NullGroup,
         typename T3=NullGroup,typename T3DB=NullGroup,
         typename T4=NullGroup,
         typename T5=NullGroup>
class AddGroupFunctor:public LibMoleculeBase{
   public:
      static void AddGroup(int Index,std::vector<bool>& IsAssigned,
            const Conn_t& Conns, SharedMol Mol,PsiMap_t& FoundGroups);
};

template<int N,typename Base,typename T1,typename T1DB, typename T1TB>
class UpdateGroupFunctor{
   public:
      static void UpdateGroup(int size,PsiMap_t& FoundGroups,
                              std::vector<int>& IDs){
         switch(size){
            case(N):{
               UpdateGroups(FoundGroups,T1(&IDs[0]));
               break;
            }
            case(N-1):{
               UpdateGroups(FoundGroups,T1DB(&IDs[0]));
               break;
            }
            case(N-2):{
               UpdateGroups(FoundGroups,T1TB(&IDs[0]));
               break;
            }
            default:{
               UpdateGroups(FoundGroups,Base(&IDs[0]));
               break;
            }
         }
      }
};

template<int N,typename Base,typename T1,typename T1DB,typename T1TB,
         typename T2,typename T2DB,typename T2TB,
         typename T3,typename T3DB,
         typename T4,
         typename T5>
void AddGroupFunctor<N,Base,T1,T1DB,T1TB,T2,T2DB,T2TB,T3,T3DB,T4,T5>::
            AddGroup(int Index,std::vector<bool>& IsAssigned,
      const Conn_t& Conns, SharedMol Mol,PsiMap_t& FoundGroups){
   //We now know index's group (or we will be by time we finish this)
   IsAssigned[Index]=true; //
   std::vector<int> IDs;
   IDs.push_back(Index);
   int NBonds=Conns[Index].size();
   typedef std::vector<int>::const_iterator cIntIt;
   for (cIntIt BondI=Conns[Index].begin(); BondI!=Conns[Index].end(); ++BondI) {
      int ConAtom=(*BondI);
      if (!IsAssigned[ConAtom]) {
         if ((*Mol)[ConAtom]->Z()==1) {
            IDs.push_back(ConAtom);
            IsAssigned[ConAtom]=true;
         }
      }
   }
   switch(IDs.size()){
      case(1):{
         UpdateGroupFunctor<N,Base,T1,T1DB,T1TB>::
             UpdateGroup(Conns[IDs[0]].size(),FoundGroups,IDs);
         break;
      }
      case(2):{
         UpdateGroupFunctor<N,Base,T2,T2DB,T2TB>::
             UpdateGroup(Conns[IDs[0]].size(),FoundGroups,IDs);
         break;
      }
      case(3):{
         UpdateGroupFunctor<N,Base,T3,T3DB,NullGroup>::
             UpdateGroup(Conns[IDs[0]].size(),FoundGroups,IDs);
         break;
      }
      case(4):{
         UpdateGroups(FoundGroups,T4(&IDs[0]));
         break;
      }
      case(5):{
         UpdateGroups(FoundGroups,T5(&IDs[0]));
         break;
      }

   }
}

}}


#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDGROUPFUNCTOR_H_ */
