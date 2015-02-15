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

#include "Fragmenter.h"
#include "LibFragFragment.h"
#include "LibFragMolecule.h"
#include "FxnalGroup.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<const Molecule> SharedMol;
typedef boost::shared_ptr<Fragment> SharedFrag;
typedef boost::shared_ptr<const FxnalGroup> SharedGroup;
typedef std::vector<SharedGroup> vSharedGroup;
BondFragmenter::BondFragmenter(SharedMol Mol, const int NBonds):
      Fragmenter(Mol),Geom_(new OrganicGeom(Mol.get())),
            NBonds_(NBonds){
}

void BondFragmenter::AddFragment(const vSharedGroup& FoundGroups,
      const long int value){
   SharedFrag temp(new Fragment(Mol_,value));
   vSharedGroup::const_iterator GroupI=FoundGroups.begin(),
                                GroupEnd=FoundGroups.end();
   for(;GroupI!=GroupEnd;++GroupI){
      for(int i=0;i<(*GroupI)->size();i++)
         (*temp)<<(*(*GroupI))[i];
   }
   //Now make sure our new frag isn't a subset of another
   std::vector<SharedFrag>::iterator FragI=FoundFrags_.begin(),
                                        FragEnd=FoundFrags_.end();
   bool good=true;
   for(;FragI!=FragEnd;++FragI){
      if((*temp)<=(*(*FragI))){
         good=false;
         break;
      }
      else if((*(*FragI))>(*temp)){
         good=false;
         (*FragI)=temp;
      }
   }
   if(good)FoundFrags_.push_back(temp);
}

void BondFragmenter::Recurse(vSharedGroup& FoundGroups,
            const Connections& Conns,
            const ConnGroups& FxnGroups,
            long int& value){
   if(FoundGroups.size()==NBonds_+1||
      FoundGroups.back()->NAttachPoint()==0){
      AddFragment(FoundGroups,value++);
      return;
   }
   //Will be true if we find a good group to follow
   bool FoundGroup=false;
   for(int i=0;i<FoundGroups.back()->NAttachPoint();i++){
      int NAttachI=FoundGroups.back()->AttachPoint(i);
      if(NAttachI==-1)continue;
      for(int j=0;j<Conns[NAttachI].size();j++){
         bool good=true;
         int ConnJ=Conns[NAttachI][j];
         //Make sure atom is an attachment point for some group
         if(FxnGroups.count(ConnJ)!=1)continue;
         //Make sure we don't have this group already
         for(int k=0;k<FoundGroups.size()-1&&good;k++){
            int NAttachK=FoundGroups[k]->NAttachPoint();
            for(int l=0;l<NAttachK&&good;l++)
               if(FoundGroups[k]->AttachPoint(l)==ConnJ)good=false;
         }
         if(good){
            FoundGroup=true;
            FoundGroups.push_back(FxnGroups.AttachedGroup(ConnJ));
            Recurse(FoundGroups,Conns,FxnGroups,value);
            FoundGroups.erase(FoundGroups.end()-1,FoundGroups.end());
         }
      }
   }
   if(!FoundGroup)//We hit a dead end
      AddFragment(FoundGroups,value++);

}

std::vector<SharedFrag> BondFragmenter::MakeFrags(){
   const ConnGroups& FxnGroups=Geom_->GetGroups();
   GroupIt FxnGroupI=FxnGroups.begin(),EndGroup=FxnGroups.end();
   const Connections& Conns=Geom_->GetConns();
   for (; FxnGroupI!=EndGroup; ++FxnGroupI){
      vSharedGroup Groups;
      Groups.push_back(*FxnGroupI);
      long int value=0;
      Recurse(Groups,Conns,FxnGroups,value);
   }
   return FoundFrags_;
}


}}//End namespaces
