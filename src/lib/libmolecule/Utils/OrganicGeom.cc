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

typedef PrimRunner<Methane, Methyl, Methene, Methyne, C4, Alkenyl1, Alkenyl2,
      Alkenyl3, Alkynyl1, Alkynyl2> FindCarbon;
typedef PrimRunner<Ammonia, Amine1, Amine2, Amine3, Azo1, Azo2, NTB, Ammonium,
      Ammonium1, Ammonium2, Ammonium3, Ammonium4> FindNitrogen;
typedef PrimRunner<Water, Hydroxyl, Ether, ODB> FindOxygen;
typedef PrimRunner<HydrogenSulfide, Thiol, Sulfide, SDB> FindSulfur;
typedef SetRunner<FindCarbon,FindOxygen,
               FindNitrogen,FindSulfur> FindPrimitiveNodes;

/*
typedef PrimRunner<Ethene, DBCC1, DBCC2, DBCC2G, DBCC3, DBCC4> FindCC2XBond;
typedef PrimRunner<Acetylene, Ethynyl, CCTB2> FindCC3XBond;
typedef PrimRunner<Formaldimine, DBCN1, Aldimine1, Aldimine2, Ketimine1,
      Ketimine2> FindCN2XBond;
typedef PrimRunner<HydrogenCyanide, Nitrile> FindCN3XBond;
typedef PrimRunner<Formaldehyde, Aldehyde, Carbonyl> FindCO2XBond;
typedef PrimRunner<Amide1C, Amide1N,Amide2, Amide2G, Amide3> FindAmide;
typedef PrimRunner<Carboxylate, Carboxyl, Carbonate, Ester, Imide1C, Imide1N, Imide2, Imide2C, Imide3> FindCOGroups;
typedef PrimRunner<Methanol, Methoxy, HydrogenPeroxide, Hydroperoxide, Peroxide,
      Hemiacetal, Hemiketal, Acetal, Ketal, Orthoester> FindEtherGroups;
typedef PrimRunner<Azide, Diazene, DBNN1, DBNN2, Cyanate, Isocyanate, Nitro,
      Nitrate, Nitroso, Nitrite> FindNGroups;
typedef PrimRunner<Phenyl, OrthoBenzene, MetaBenzene, ParaBenzene> FindComplexRing;
typedef PrimRunner<Benzene, Benzene6> Find6Ring;
typedef PrimRunner<Isopropyl> FindMiscCGroups;

typedef PrimRunner<Indolyl0, Indolyl1_1, Indolyl1_2, Indolyl1_3, Indolyl1_5,
      Indolyl1_6, Indolyl1_7, Indolyl1_8, Indolyl2_12, Indolyl2_13, Indolyl2_15,
      Indolyl2_16, Indolyl2_17, Indolyl2_18, Indolyl2_23, Indolyl2_25,
      Indolyl2_26, Indolyl2_27, Indolyl2_28, Indolyl2_35, Indolyl2_36,
      Indolyl2_37, Indolyl2_38, Indolyl2_56, Indolyl2_57, Indolyl2_58,
      Indolyl2_67, Indolyl2_68, Indolyl2_78, Indolyl3_123, Indolyl3_125,
      Indolyl3_126, Indolyl3_127, Indolyl3_128, Indolyl3_135, Indolyl3_136,
      Indolyl3_137, Indolyl3_138, Indolyl3_156, Indolyl3_157, Indolyl3_158,
      Indolyl3_167, Indolyl3_168, Indolyl3_178, Indolyl3_235, Indolyl3_236,
      Indolyl3_237, Indolyl3_238, Indolyl3_256, Indolyl3_257, Indolyl3_258,
      Indolyl3_267, Indolyl3_268, Indolyl3_278, Indolyl3_356, Indolyl3_357,
      Indolyl3_358, Indolyl3_367, Indolyl3_368, Indolyl3_378, Indolyl3_567,
      Indolyl3_568, Indolyl3_578, Indolyl3_678, Indolyl4_1235, Indolyl4_1236,
      Indolyl4_1237, Indolyl4_1238, Indolyl4_1256, Indolyl4_1257, Indolyl4_1258,
      Indolyl4_1267, Indolyl4_1268, Indolyl4_1278, Indolyl4_1356, Indolyl4_1357,
      Indolyl4_1358, Indolyl4_1367, Indolyl4_1368, Indolyl4_1378, Indolyl4_1567,
      Indolyl4_1568, Indolyl4_1578, Indolyl4_1678, Indolyl4_2356, Indolyl4_2357,
      Indolyl4_2358, Indolyl4_2367, Indolyl4_2368, Indolyl4_2378, Indolyl4_2567,
      Indolyl4_2568, Indolyl4_2578, Indolyl4_2678, Indolyl4_3567, Indolyl4_3568,
      Indolyl4_3578, Indolyl4_3678, Indolyl4_5678, Indolyl5_12356,
      Indolyl5_12357, Indolyl5_12358, Indolyl5_12367, Indolyl5_12368,
      Indolyl5_12378, Indolyl5_12567, Indolyl5_12568, Indolyl5_12578,
      Indolyl5_12678, Indolyl5_13567, Indolyl5_13568, Indolyl5_13578,
      Indolyl5_13678, Indolyl5_15678, Indolyl5_23567, Indolyl5_23568,
      Indolyl5_23578, Indolyl5_23678, Indolyl5_25678, Indolyl5_35678,
      Indolyl6_123567, Indolyl6_123568, Indolyl6_123578, Indolyl6_123678,
      Indolyl6_125678, Indolyl6_135678, Indolyl6_235678, Indolyl7> FindIndole;

typedef SetRunner<FindIndole,FindComplexRing,FindMiscCGroups,FindCN2XBond,
      FindCN3XBond,FindCO2XBond,FindEtherGroups,FindCOGroups,FindNGroups,
      FindCC2XBond,FindCC3XBond> FindDerivedNodes;*/

Graph OrganicGeom::MakeFxnGroups(bool FindAAs) const {
   MolItr AtomI=Mol_->Begin(),AtomEnd=Mol_->End();
   std::vector<bool> IsAssigned(Mol_->NAtoms(), false);
   Graph Nodes;
   PeriodicTable PTable;
   std::vector<boost::shared_ptr<Node> > Temp;
   for (int counter=0; AtomI!=AtomEnd; ++AtomI) {
      PTable(counter++, (*AtomI), Nodes);
      Temp.push_back(Nodes.back());
   }
   //Don't want the references hanging around
   Temp.clear();
   Graph::iterator It=Nodes.begin(),ItEnd=Nodes.end();
   for (unsigned counter=0; It!=ItEnd; ++It, ++counter) {
      for (unsigned ConnI=0; ConnI<Connections_[counter].size(); ConnI++)
         (*It)->AddConn(Temp[Connections_[counter][ConnI]]);
   }
   FindPrimitiveNodes::Run(Nodes);
   SetRunner<PrimRunner<DBCC,TBCC> >::Run(Nodes);
   /*FindDerivedNodes::Run(Nodes);
   //If this is true we will find the peptide bond as an amide
   if(!FindAAs)SetRunner<FindAmide>::Run(Nodes);
   SetRunner<Find6Ring>::Run(Nodes);*/

   return Nodes;
}

OrganicGeom::OrganicGeom(const Molecule* Mol,bool FindAAs) :
      Geometry(Mol) {
   FxnalGroups_=MakeFxnGroups(FindAAs);
}

std::string OrganicGeom::PrintOut() const {
   return FxnalGroups_.PrintOut();
}
}} //End namespaces
