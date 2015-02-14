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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDDERIVEDGRPFUNCTOR_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDDERIVEDGRPFUNCTOR_H_
#include <stack>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "FxnalGroup.h"
#include "../../libPsiUtil/PsiMap.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<FxnalGroup> SharedGroup;
typedef boost::shared_ptr<const FxnalGroup> cSharedGroup;
typedef std::vector<std::vector<int> > Conn_t;
typedef std::vector<int>::const_iterator cIntIt;

template<typename T>
class AddDerivedGrp{
   public:
      static bool AddGroup(std::stack<cSharedGroup>& DaGroups,
                            const Conn_t& Conns,
                            PsiMap<int,cSharedGroup>& FoundGroups);

};

template<>
class AddDerivedGrp<NullGroup>{
   public:
      static bool AddGroup(std::stack<cSharedGroup>&,
                            const Conn_t&,
                            PsiMap<int,cSharedGroup>&){
         return false;
      }
};

template<typename T>
bool AddDerivedGrp<T>::AddGroup(std::stack<cSharedGroup>& Groups,
             const Conn_t& Conns,
             PsiMap<int,cSharedGroup>& FoundGroups){
   bool found=false;
   T TempObj;
   if(Groups.size()!=TempObj.GetNPrims()){
      FxnGroup_t GroupType=TempObj.GetPrimI(Groups.size());
      int AtomWeWant=Groups.top()->AttachPoint();
      cIntIt AtomI=Conns[AtomWeWant].begin(),
            LastAtom=Conns[AtomWeWant].end();
      for(;AtomI!=LastAtom;++AtomI){
         if(FoundGroups.count(*AtomI)==1){
            cSharedGroup AGroup=FoundGroups[(*AtomI)];
            if(AGroup->Type()==GroupType){//Have desired group
               Groups.push(AGroup);
               found=
                AddDerivedGrp<T>::AddGroup(Groups,Conns,FoundGroups);
            }
         }
      }
   }
   else{
      //Right now no derived group contains a derived group w/ more than
      //one attachment point, so we can get away only erasing the first
      //attachment point
       std::vector<cSharedGroup> NewGroup;
       while(!Groups.empty()){
          cSharedGroup Last=Groups.top();
          Groups.pop();
          FoundGroups.erase(Last->AttachPoint());
          NewGroup.push_back(Last);
       }
       boost::shared_ptr<T> temp(new T(NewGroup));
       FoundGroups[NewGroup[0]->AttachPoint(0)]=temp;
       for(int i=1;i<temp->NAttachPoint();i++)
          FoundGroups[temp->AttachPoint(i)]=temp;
       found=true;
   }
   return found;
}

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ADDDERIVEDGRPFUNCTOR_H_ */
