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
#include "OrganicGeom.h"
#include "BioGeom.h"
#include "AutoFxnalGroup.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<const Molecule> SharedMol;
typedef boost::shared_ptr<Fragment> SharedFrag;
typedef boost::shared_ptr<const Node> SharedGroup;
typedef std::vector<SharedGroup> vSharedGroup;
BondFragmenter::BondFragmenter(SharedMol Mol, const unsigned int NBonds):
      Fragmenter(Mol),Geom_(new OrganicGeom(Mol.get())),
            NBonds_(NBonds){
   //BioGeom BChem(*Mol);
}

void BondFragmenter::AddFragment(const vSharedGroup& FoundGroups,
      const long int value){
   SerialNumber TempSN;
   TempSN.insert(value+1);
   SharedFrag temp(new Fragment(Mol_,TempSN));
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
            const Graph& FxnGroups,
            long int& value){
   if(FoundGroups.size()==NBonds_+1
      //FoundGroups.back()->NAttachPoint()==0
         ){
      AddFragment(FoundGroups,value++);
      return;
   }
   //Will be true if we find a good group to follow
   bool FoundGroup=false;
   for(int i=0;i<FoundGroups.back()->ConnNodes_.size();i++){
      SharedGroup NodeI=FoundGroups.back()->ConnNodes_[i];
      bool good=true;
      //Make sure we don't have this group already
      for(unsigned k=0;k<FoundGroups.size()-1&&good;k++)
         if(FoundGroups[k].get()==NodeI.get())break;
      if(good){
         FoundGroup=true;
         FoundGroups.push_back(NodeI);
         Recurse(FoundGroups,Conns,FxnGroups,value);
         FoundGroups.erase(FoundGroups.end()-1,FoundGroups.end());
      }
   }
   if(!FoundGroup)//We hit a dead end
      AddFragment(FoundGroups,value++);
}

std::vector<SharedFrag> BondFragmenter::MakeFrags(){
   const Graph& FxnGroups=Geom_->GetGroups();
   std::cout<<FxnGroups.PrintOut();
   Graph::const_iterator FxnGroupI=FxnGroups.begin(),EndGroup=FxnGroups.end();
   const Connections& Conns=Geom_->GetConns();
   long int value=0;
   for (; FxnGroupI!=EndGroup; ++FxnGroupI){
      vSharedGroup Groups;
      Groups.push_back(*FxnGroupI);
      Recurse(Groups,Conns,FxnGroups,value);
   }
   return FoundFrags_;
}


}}//End namespaces
