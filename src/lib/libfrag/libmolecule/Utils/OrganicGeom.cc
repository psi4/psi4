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
#include "AssignAtomFunctor.h"
#include "FindDerGroup.h"
#include "AddGroupFunctor.h"
#include "AddDerivedGrpFunctor.h"
#include "AromaRingFinder.h"

namespace psi{
namespace LibMolecule{

template<typename T,typename...Args>
class PrimRunner: public PrimRunner<Args...>{
   public:
      static void Run(GroupIt& Group,const Connections& Conns,
                      ConnGroups& FoundGroups){
         if(T::FindMe::FindGroup((*Group),Conns,FoundGroups))
            Group=FoundGroups.begin();
         else
            PrimRunner<Args...>::Run(Group,Conns,FoundGroups);
      }
};

template<typename T>
class PrimRunner<T>{
   public:
      static void Run(GroupIt& Group,const Connections& Conns,
            ConnGroups& FoundGroups){
         if(T::FindMe::FindGroup((*Group),Conns,FoundGroups))
            Group=FoundGroups.begin();
      }
};

typedef PrimRunner<Water,Hydroxyl,Oxygen2,OxyDB> OxygenPrims;
typedef PrimRunner<Methane,Methyl,Methene,Methyne,Carbon4,
                   Alkenyl1,Alkenyl2,Alkenyl3,
                   Alkynyl1,Alkynyl2> CarbonPrims;
typedef PrimRunner<Ammonia,Amine1,Amine2,Amine3,
                   NitroDB1,NitroDB2,NitroTB> NitrogenPrims;
typedef PrimRunner<HydrogenFluoride,Fluorine1> FlourinePrims;
typedef PrimRunner<HydrogenChloride,Chlorine1> ChlorinePrims;
typedef PrimRunner<HydrogenBromide,Bromine1> BrominePrims;
typedef PrimRunner<HydrogenIodide,Iodine1> IodinePrims;

typedef PrimRunner<HydrogenCyanide,Ethyne> Alkynyl1Groups;
typedef PrimRunner<Nitrile,cctb,Ethynyl>Alkynyl2Groups;
typedef PrimRunner<Formaldehyde,Ethene,Methanimine> Alkenyl1Groups;
typedef PrimRunner<Aldehyde,ccdb3,ccdb2,Ethenyl1,Aldimine1,Aldimine2> Alkenyl2Groups;
typedef PrimRunner<Carbonyl,ccdb4,Ethenyl2,Ketimine1,Ketimine2> Alkenyl3Groups;
typedef PrimRunner<Carboxyl> CarbonylGroups;
typedef PrimRunner<Methoxy,HydroPeroxy> Oxygen2Groups;
typedef PrimRunner<Peroxide,Methanol> HydroxylGroups;

ConnGroups  OrganicGeom::MakeFxnGroups()const {
   MolItr AtomI=Mol_->Begin(),AtomEnd=Mol_->End();
   std::vector<bool> IsAssigned(Mol_->NAtoms(), false);
   ConnGroups FoundGroups;
   PeriodicTable PTable;
   for(int counter=0;AtomI!=AtomEnd;++AtomI)
      PTable(counter++,(*AtomI),FoundGroups);
   GroupIt GroupI=FoundGroups.begin();
   for(;GroupI!=FoundGroups.end();++GroupI){
      if((*GroupI)->Type()==CARBON)
               CarbonPrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==NITROGEN)
         NitrogenPrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==OXYGEN)
         OxygenPrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==FLUORINE)
         FlourinePrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==CHLORINE)
         ChlorinePrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==BROMINE)
         BrominePrims::Run(GroupI,Connections_,FoundGroups);
      else if((*GroupI)->Type()==IODINE)
         IodinePrims::Run(GroupI,Connections_,FoundGroups);
   }
   for(GroupI=FoundGroups.begin();GroupI!=FoundGroups.end();++GroupI){
	   Alkynyl1Groups::Run(GroupI,Connections_,FoundGroups);
	   Alkynyl2Groups::Run(GroupI,Connections_,FoundGroups);
	   Alkenyl1Groups::Run(GroupI,Connections_,FoundGroups);
	   Alkenyl2Groups::Run(GroupI,Connections_,FoundGroups);
	   Alkenyl3Groups::Run(GroupI,Connections_,FoundGroups);
	   CarbonylGroups::Run(GroupI,Connections_,FoundGroups);
	   Oxygen2Groups::Run(GroupI,Connections_,FoundGroups);
	   HydroxylGroups::Run(GroupI,Connections_,FoundGroups);
   }
   std::cout<<FoundGroups.PrintOut();
   exit(1);
   /*/Step 1: Find basic groups
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
      found=Finder.FindRing(FoundGroups,Connections_);*/
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

std::string ConnGroups::PrintOut()const{
   GroupIt GroupI=begin(),GroupEnd=end();
   std::stringstream Result;
   for(;GroupI!=GroupEnd;++GroupI)
      Result<<(*GroupI)->PrintOut("");
   return Result.str();
}

std::string OrganicGeom::PrintOut()const{
   return  FxnalGroups_.PrintOut();
}
}}//End namespaces
