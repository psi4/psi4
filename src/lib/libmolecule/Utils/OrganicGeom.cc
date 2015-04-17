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
#include "AutoFxnalGroup.h"
#include "AutoFxnalGroup/PrimRunner.h"
namespace psi {
namespace LibMolecule {

///A simple functor to assign atoms to fxnalgroups based on Z
template <int Z, typename T, typename ...Args>
class AssignAtomFunctor:public AssignAtomFunctor<Z+1, Args...> {
   private:
      typedef AssignAtomFunctor<Z+1, Args...> Base_t;
   public:
      void operator()(int i, const boost::shared_ptr<const Atom> AtomIn,
            Graph& Groups) const {
         if (AtomIn->Z()==Z) Groups.push_back(boost::shared_ptr<T>(new T(i)));
         else Base_t::operator()(i, AtomIn, Groups);
      }
};

///Termination condition to previous functor, prevents infinite recursion
template <int Z, typename T>
class AssignAtomFunctor<Z, T> {
   public:
      void operator()(int i, const boost::shared_ptr<const Atom> AtomIn,
            Graph& Groups) const {
         if (AtomIn->Z()==Z) Groups.push_back(boost::shared_ptr<T>(new T(i)));

      }
};

typedef AssignAtomFunctor<1,Hydrogen,Helium,Lithium,Beryllium,
   Boron,Carbon,Nitrogen,Oxygen,Fluorine,Neon,Sodium,Magnesium,Aluminum,
   Silicon,Phosphorus,Sulfur,Chlorine,Argon,Potassium,Calcium,Scandium,
   Titanium,Vanadium,Chromium,Manganese,Iron,Cobalt,Nickel,Copper,Zinc,
   Gallium,Germanium,Arsenic,Selenium,Bromine,Krypton,Rubidium,Strontium,
   Yttrium,Zirconium,Niobium,Molybdenum,Technetium,Ruthenium,Rhodium,
   Palladium,Silver,Cadmium,Indium,Tin,Antimony,Tellurium,Iodine,Xenon,
   Cesium,Barium,Lanthanum,Cerium,Praseodymium,Neodymium,Promethium,
   Samarium,Europium,Gadolinium,Terbium,Dysprosium,Holmium,Erbium,Thulium,
   Ytterbium,Lutetium,Hafnium,Tantalum,Tungsten,Rhenium,Osmium,Iridium,
   Platinum,Gold,Mercury,Thallium,Lead,Bismuth,Polonium,Astatine,Radon,
   Francium,Radium,Actinium,Thorium,Protactinium,Uranium,Neptunium,
   Plutonium,Americium,Curium,Berkelium,Californium,Einsteinium,Fermium,
   Mendelevium,Nobelium,Lawrencium,Rutherfordium,Dubnium,Seaborgium,Bohrium,
   Hassium,Meitnerium,Darmstadtium,Roentgenium,Copernicium,Ununtrium,
   Flerovium,Ununpentium,Livermorium,Ununseptium,Ununoctium> PeriodicTable;

typedef PrimRunner<Methane, Methyl, Methene,
                   Methyne, C4, Alkenyl1, Alkenyl2,
                   Alkenyl3, Alkynyl1, Alkynyl2> FindCarbon;
typedef PrimRunner<Ammonia, Amine1, Amine2, Amine3,
                   Azo1, Azo2, NTB, Ammonium,
                   Ammonium1, Ammonium2,
                   Ammonium3, Ammonium4> FindNitrogen;
typedef PrimRunner<Water, Hydroxyl, Ether, ODB> FindOxygen;
typedef PrimRunner<HydrogenSulfide, Thiol, Sulfide, SDB> FindSulfur;
typedef SetRunner<FindCarbon,FindOxygen,
               FindNitrogen,FindSulfur> FindPrimitiveNodes;

typedef PrimRunner<Ethene, DBCC> FindCC2XBond;
typedef PrimRunner<Acetylene, TBCC> FindCC3XBond;
typedef PrimRunner<Formaldimine, DBCN1, Aldimine1, Aldimine2, Ketimine1,
      Ketimine2> FindCN2XBond;
typedef PrimRunner<HydrogenCyanide, Nitrile> FindCN3XBond;
typedef PrimRunner<Formaldehyde, Aldehyde, Carbonyl> FindCO2XBond;
typedef PrimRunner<Imide,Amide> FindAmide;
typedef PrimRunner<Carboxylate, Carboxyl, Carbonate, Ester> FindCOGroups;
typedef PrimRunner<Methanol, Methoxy, HydrogenPeroxide, Hydroperoxide, Peroxide,
       OrthoEster> FindEtherGroups;
typedef PrimRunner<Azide, Diazene, DBNN, Cyanate, Isocyanate, Nitro,
      Nitrate, Nitroso, Nitrite> FindNGroups;
typedef PrimRunner<Hemiacetal,Hemiketal,Acetal,Ketal> FindAcetal;
typedef SetRunner<FindCN2XBond,
      FindCN3XBond,FindCO2XBond,FindEtherGroups,FindCOGroups,FindNGroups,
      FindCC2XBond,FindCC3XBond> FindDerivedNodes;

typedef PrimRunner<Benzofuran,Isobenzofuran,Indole,
      Isoindole,Purine,Indazole,Benzoxazole,Benzisoxazole> Find9MemberRing;
typedef PrimRunner<Furan,Pyrrole,Thiophene,Imidazole,
      Pyrazole,Oxazole,Isoxazole> Find5MemberRing;
typedef PrimRunner<Benzene,Pyridine,Pyrazine,Pyrimidine,Pyridazine>
      Find6MemberRing;
typedef SetRunner<Find6MemberRing,Find5MemberRing> FindRingNodes;

Graph OrganicGeom::MakeFxnGroups(bool FindAAs) const {
   MolItr AtomI=Mol_.Begin(),AtomEnd=Mol_.End();
   std::vector<bool> IsAssigned(Mol_.NAtoms(), false);
   Graph Nodes;
   PeriodicTable PTable;
   std::vector<boost::shared_ptr<Node> > Temp;
   for (int counter=0; AtomI!=AtomEnd; ++AtomI) {
      PTable(counter++, (*AtomI), Nodes);
      Temp.push_back(Nodes.back());
   }
   Graph::iterator It=Nodes.begin(),ItEnd=Nodes.end();
   for (unsigned counter=0; It!=ItEnd; ++It, ++counter) {
      for (unsigned ConnI=0; ConnI<Connections_[counter].size(); ConnI++)
         (*It)->AddConn(Temp[Connections_[counter][ConnI]],*It);
   }
   //Don't want the references hanging around
   Temp.clear();
   FindPrimitiveNodes::Run(Nodes);
   SetRunner<Find9MemberRing>::Run(Nodes);
   FindRingNodes::Run(Nodes);
   FindDerivedNodes::Run(Nodes);
   if(!FindAAs)SetRunner<FindAmide>::Run(Nodes);
   //Acetal groups really get in the way of other groups so we look for them
   //last
   SetRunner<FindAcetal>::Run(Nodes);
   /*If this is true we will find the peptide bond as an amide
   if(!FindAAs)SetRunner<FindAmide>::Run(Nodes);*/

   return Nodes;
}

OrganicGeom::OrganicGeom(const Molecule& Mol,bool FindAAs) :
      Geometry(Mol) {
   FxnalGroups_=MakeFxnGroups(FindAAs);
}

std::string OrganicGeom::PrintOut() const {
   return FxnalGroups_.PrintOut();
}
}} //End namespaces
