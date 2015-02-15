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
#include "OrganicGeom.h"
#include "MolItr.h"
#include "LibFragMolecule.h"
#include "FxnalGroup.h"
#include "AddGroupFunctor.h"
#include "AddDerivedGrpFunctor.h"
#include "AromaRingFinder.h"
namespace psi{
namespace LibMolecule{
typedef AddGroupFunctor<4,Carbon,Carbon4,Alkenyl3,Alkynyl2,
                        Methyne,Alkenyl2,Alkynyl1,
                        Methene,Alkenyl1,
                        Methyl,
                        Methane> CFunctor;
typedef AddGroupFunctor<2,Oxygen,Oxygen2,OxyDB,NullGroup,
                       Hydroxyl,NullGroup,NullGroup,
                       Water> OFunctor;
typedef AddGroupFunctor<3,Nitrogen,Amine3,NitroDB2,NitroTB,
                        Amine2,NitroDB1,NullGroup,
                        Amine1,NullGroup,
                        Ammonia> NFunctor;
typedef AddGroupFunctor<1,Fluorine,Fluorine1,NullGroup,NullGroup,
                        HydrogenFluoride> FFunctor;
typedef AddGroupFunctor<1,Chlorine,Chlorine1,NullGroup,NullGroup,
                        HydrogenChloride> ClFunctor;
typedef AddGroupFunctor<1,Bromine,Bromine1,NullGroup,NullGroup,
                        HydrogenBromide>   BrFunctor;
typedef AddGroupFunctor<1,Iodine,Iodine1,NullGroup,NullGroup,
                        HydrogenIodide>     IFunctor;

template<typename T1,typename T2=NullGroup,
         typename T3=NullGroup,typename T4=NullGroup>
class DerGrpFunctor{
   public:
      static bool AddGroup(
            std::stack<boost::shared_ptr<const FxnalGroup> >& DaGroups,
            const std::vector<std::vector<int> >& Connections,
            PsiMap_t& FoundGroups){
         bool LoopRestart=AddDerivedGrp<T1>::AddGroup(
               DaGroups,Connections,FoundGroups);
         if(!LoopRestart)
            LoopRestart=AddDerivedGrp<T2>::AddGroup(
                  DaGroups,Connections,FoundGroups);
         if(!LoopRestart)
            LoopRestart=AddDerivedGrp<T3>::AddGroup(
                  DaGroups,Connections,FoundGroups);
         if(!LoopRestart)
            LoopRestart=AddDerivedGrp<T4>::AddGroup(
                  DaGroups,Connections,FoundGroups);
         return LoopRestart;
      }
};


ConnGroups  OrganicGeom::MakeFxnGroups()const {
   MolItr AtomI=Mol_->Begin();
   std::vector<bool> IsAssigned(Mol_->NAtoms(), false);
   ConnGroups FoundGroups;
   //Step 1: Find basic groups
   for (int index=0; AtomI!=Mol_->End(); ++AtomI, ++index) {
      if (IsAssigned[index]) continue;
      switch (AtomI->Z()) {
         case (6): {
            CFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case (7): {
            NFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case (8): {
            OFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case(9):{
            FFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case(17):{
            ClFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case(35):{
            BrFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }
         case(53):{
            IFunctor::AddGroup(index,IsAssigned,Connections_,Mol_,FoundGroups);
            break;
         }

      }

   }
   //Step 2: Find derived groups
   bool AllFound=false;
   while (!AllFound) {
      bool LoopRestart=false;
      GroupIt FGrpI=FoundGroups.begin();
      for (; FGrpI!=FoundGroups.end(); ++FGrpI) {
         std::stack<boost::shared_ptr<const FxnalGroup> > DaGroups;
         DaGroups.push((*FGrpI));
         switch ((*FGrpI)->Type()) {
            case (NITROTB): {
               LoopRestart=
                 DerGrpFunctor<Nitrile,HydrogenCyanide>::
                    AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case (OXYDB):{
               //These groups all assume O is double bonded to something
               LoopRestart=
                 DerGrpFunctor<Aldehyde,Carbonyl,Formaldehyde>::
                    AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case (HYDROXYL):{
               LoopRestart=
                 DerGrpFunctor<Carboxyl,HydroPeroxy,Methanol,Peroxide>::
                    AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case (ALKENYL1):{
               LoopRestart=
                     DerGrpFunctor<Ethenyl2,Ethenyl1,Ethene>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case(ALKENYL2):{
               LoopRestart=DerGrpFunctor<ccdb3,ccdb2>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case(ALKENYL3):{
               LoopRestart=DerGrpFunctor<ccdb4>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case (ALKYNYL1):{
               LoopRestart=
                     DerGrpFunctor<Ethynyl,Ethyne>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case(ALKYNYL2):{
               LoopRestart=DerGrpFunctor<cctb>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case(NITRODB1):{
               LoopRestart=DerGrpFunctor<Aldimine1,Ketimine1,Methanimine>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
            case(NITRODB2):{
               LoopRestart=DerGrpFunctor<Aldimine2,Ketimine2>::
                     AddGroup(DaGroups,Connections_,FoundGroups);
               break;
            }
         }
         if(LoopRestart)break;
      }
      if (!LoopRestart) AllFound=true;
   }
   AromaticRingFinder Finder;
   bool found=true;
   while(found)
      found=Finder.FindRing(FoundGroups,Connections_);
   return FoundGroups;
}

OrganicGeom::OrganicGeom(const Molecule* Mol,const double BondDef,
      const int MaxBonds):
   Geometry(Mol,BondDef,MaxBonds){
   FxnalGroups_=MakeFxnGroups();
}

GroupIt& GroupIt::operator++(){
   //Increment is outside loop to ensure it happens at least once
   ++It_;
   for(;It_!=ItEnd_;++It_)
      //If our 1st increment was good, we never hit the loop increment
      if((It_->first)==((It_->second)->AttachPoint()))break;
   return *this;
}

std::string OrganicGeom::PrintOut()const{
   GroupIt GroupI=FxnalGroups_.begin(),GroupEnd=FxnalGroups_.end();
   std::stringstream Result;
   for(;GroupI!=GroupEnd;++GroupI)
      Result<<(*GroupI)->PrintOut();
   return Result.str();
}
}}//End namespaces
