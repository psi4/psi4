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
#include "../libpsio/MOFile.h"
#include "../libmints/matrix.h"
#include "../libmints/basisset_parser.h"
#include "../libmints/basisset.h"
#include "MBEFragSet.h"
#include "MBEFrag.h"
#include <vector>
#include <utility>
#include <map>
namespace psi {
namespace LibFrag {

typedef boost::shared_ptr<const MBEFrag> ShareFrag;
typedef std::map<std::string, SharedMatrix> SpinMats;
typedef std::map<std::string, int> SpinDims;
typedef std::vector<int> IntVec;
typedef std::map<int, int> IntMap;


inline void DimsPerAtom(IntVec& Atoms, IntMap& AtomRowStart,
      IntMap& AtomColStart, IntMap& ElecPAtom) {
   boost::shared_ptr<BasisSet> primary=BasisSet::pyconstruct_orbital(
         Process::environment.molecule(),
         "BASIS", Process::environment.options.get_str("BASIS"));
   boost::shared_ptr<Molecule> Mol=primary->molecule();
   int NAtoms=Atoms.size();
   for (int i=0,RowCount=0,ColCount=0; i<NAtoms; i++) {
      int CurrAtom=Atoms[i];
      AtomRowStart[CurrAtom]=RowCount;
      AtomColStart[CurrAtom]=ColCount;
      for (int shell=0; shell<primary->nshell_on_center(CurrAtom); shell++) {
         const GaussianShell& Shell=primary->shell(CurrAtom, shell);
         RowCount+=Shell.nfunction();
      }
      ElecPAtom[CurrAtom]=Mol->Z(CurrAtom);
      ColCount+=Mol->Z(CurrAtom);
   }
}

inline MOFile WriteFile(double energy, int NBf, SpinDims& ColDims,
      SpinMats& Cs) {
   boost::shared_ptr<BasisSet> primary=BasisSet::pyconstruct_orbital(
         Process::environment.molecule(),
         "BASIS", Process::environment.options.get_str("BASIS"));
   MOFile NewFile;
   int length=primary->name().length()+1;
   NewFile.FillFile(length, primary->name().c_str(), primary->has_puream(),
         energy, 1, NBf, &NBf, &ColDims["ALPHA"], &ColDims["BETA"], Cs["ALPHA"],
         Cs["BETA"]);
   return NewFile;
}

MOFile DirectSum(const int N, const int x, const std::vector<MOFile>& MOFiles,
      const std::vector<MBEFragSet>& Systems) {
   std::vector<std::string>Spins;
   Spins.push_back("ALPHA");
   Spins.push_back("BETA");
   //These will be the total dimensions
   SpinDims RowDims,ColDims;
   for (int i=0; i<Spins.size(); i++) {
      RowDims[Spins[i]]=0;
      ColDims[Spins[i]]=0;
   }

   //These are the alpha and beta Cs by frag
   std::vector<SpinMats> Coefs;

   //Our Nmer
   ShareFrag NMer=Systems[N][x];

   ///This will hold all atoms in the NMer and in the final order (sorted)
   IntVec Atoms;

   ///A map between overall atom number and overall(absolute) fragment
   IntMap Atom2Frag;

   ///Starting point of each fragment's orbitals by spin, indexed by absolute
   ///frag number
   std::map<int,SpinDims> OrbStartsPFrag;

   ///I need to be able to convert one to the other b/c MOFiles are absolute
   std::map<int,int> RelFrag2AbsFrag;
   std::map<int,int> AbsFrag2RelFrag;

   ///The sum of all the energies
   double energy=0;

   for (int frag=0; frag<N+1; frag++) {
      int Fragi=NMer->ParentI(frag);
      RelFrag2AbsFrag[frag]=Fragi;
      AbsFrag2RelFrag[Fragi]=frag;
      ShareFrag MonoI=Systems[0][Fragi];
      for (int j=0; j<MonoI->Atoms().size(); j++) {
         Atoms.push_back(MonoI->Atoms()[j]);
         Atom2Frag[Atoms.back()]=Fragi;
      }
      Coefs.push_back(SpinMats());
      OrbStartsPFrag[Fragi]=SpinDims();
      energy+=MOFiles[Fragi].GetEnergy();
      for (int spin=0; spin<2; spin++) {
         Coefs.back()[Spins[spin]]=(
               Spins[spin]=="ALPHA" ? MOFiles[Fragi].GetCa() :
                                      MOFiles[Fragi].GetCb());
         RowDims[Spins[spin]]+=Coefs.back()[Spins[spin]]->nrow();
         int nelecs=Coefs.back()[Spins[spin]]->ncol();
         OrbStartsPFrag[Fragi][Spins[spin]]=ColDims[Spins[spin]];
         ColDims[Spins[spin]]+=nelecs;
      }
   }
   SpinMats Results;
   for (int i=0; i<2; i++) {
      std::string label=Spins[i]+" MOS";
      Results[Spins[i]]=SharedMatrix(
            new Matrix(label.c_str(), 1, &RowDims[Spins[i]], &ColDims[Spins[i]]));
      Results[Spins[i]]->set(0.0);
   }

   std::sort(Atoms.begin(), Atoms.end());

   //Rows are basis functions, Cols are Occupied orbitals
   IntMap AtomRowStart,AtomColStart,ElecPAtom;
   DimsPerAtom(Atoms, AtomRowStart, AtomColStart, ElecPAtom);

   std::vector<int> FragRowIndices(N+1, 0);
   int NAtoms=Atoms.size();
   for (int i=0; i<NAtoms; i++) {
      int CurrRowAtom=Atoms[i];
      int AbsFrag=Atom2Frag[CurrRowAtom];
      int RelFrag=AbsFrag2RelFrag[AbsFrag];
      //This is the number of rows our atom spans (independent of spin)
      int nrows=(i<NAtoms-1 ? AtomRowStart[Atoms[i+1]] : RowDims["ALPHA"])
                        -AtomRowStart[CurrRowAtom];
      for (int j=0; j<nrows; j++) {
         for (int spin=0; spin<2; spin++) {
            std::string SpinI=Spins[spin];
            int FullRow=AtomRowStart[CurrRowAtom]+j;
            int FragRow=FragRowIndices[RelFrag];
            int ncols=(RelFrag<N?
                     OrbStartsPFrag[RelFrag2AbsFrag[RelFrag+1]][SpinI]:
                     ColDims[SpinI])-OrbStartsPFrag[AbsFrag][SpinI];
            for (int k=0; k<ncols; k++) {
               int FullCol=OrbStartsPFrag[AbsFrag][SpinI]+k;
               int FragCol=k;
               double value=(*Coefs[RelFrag][SpinI])(FragRow, FragCol);
               Results[SpinI]->set(FullRow, FullCol, value);
            }
         }
         FragRowIndices[RelFrag]++;
      }
   }

   return WriteFile(energy, RowDims["ALPHA"], ColDims, Results);

}

}
}

