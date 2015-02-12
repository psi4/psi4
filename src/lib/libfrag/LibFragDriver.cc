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
#include <vector>
#include <boost/python/object.hpp>
#include "psi4-dec.h"
#include "process.h"
#include "LibFragDriver.h"
#include "libmolecule/LibMolecule.h"
#include "../libmints/molecule.h"
#include "libparallel2/LibParallel2.h"
#include "libPsiUtil/PythonFxn.h"
#include "libPsiUtil/Timer.h"
#include "MBEProp.h"
#include "MBE.h"
typedef boost::python::str PyStr;
namespace psi {
namespace LibFrag {
typedef boost::shared_ptr<Molecule> SharedMol;
typedef boost::shared_ptr<LibMolecule::FragmentedSystem> SharedFrags;
typedef boost::shared_ptr<LibMolecule::Fragment> SharedFrag;
typedef LibMolecule::FragmentedSystem::iterator FragItr;

SharedFrags MakeMolecule(const int N, SharedMol AMol) {
   boost::shared_ptr<LibMolecule::Molecule> MyMol(new LibMolecule::Molecule);
   for (int i=0; i<AMol->natom(); ++i) {
      std::vector<double> Carts(3, 0.0);
      Carts[0]=AMol->x(i);
      Carts[1]=AMol->y(i);
      Carts[2]=AMol->z(i);
      (*MyMol)<<LibMolecule::Atom(&Carts[0], (int)AMol->Z(i));
   }
   std::vector<double> Angles(3, 0.0),Sides(3, 0.0);
   Options options=psi::Process::environment.options;
   int SCSides=options["SUPER_CELL_SIDES"].size();
   SharedFrags MySys;
   if (SCSides>=1) {
      std::vector<int> SCellSides(3, 0);
      for (int i=0; i<3; i++) {
         Sides[i]=options["UNIT_CELL_SIDES"][i].to_double();
         Angles[i]=options["UNIT_CELL_ANGLES"][i].to_double();
         SCellSides[i]=options["SUPER_CELL_SIDES"][(SCSides==3 ? i : 0)]
               .to_integer();
      }
      boost::shared_ptr<LibMolecule::UnitCell> UC(
            new LibMolecule::UnitCell(*MyMol, &Sides[0], &Angles[0]));
      UC->FixUnitCell();
      if (SCellSides[0]+SCellSides[1]+SCellSides[2]>1) {
         boost::shared_ptr<LibMolecule::SuperCell> temp(
               new LibMolecule::SuperCell((*UC),SCellSides));
         std::cout<<temp->PrintOut(0);
         MySys=SharedFrags(
               new LibMolecule::FragmentedSystem(*temp, N));
      }
      else MyMol=UC;
   }
   if (!MySys)
      MySys=SharedFrags(new LibMolecule::FragmentedSystem(*MyMol, N));
   return MySys;
}

LibFragDriver::LibFragDriver(){
   Options options=psi::Process::environment.options;
   int N=options["TRUNCATION_ORDER"].to_integer();
   SharedMol AMol=psi::Process::environment.molecule();
   Frags_=MakeMolecule(N,AMol);
   //(*outfile)<<Frags_->PrintOut()<<std::endl;
   Energies_=boost::shared_ptr<MBEProp<double> >(
         new MBEProp<double>(N,"Energies"));
   SmartTimer MonoTimer("Run Monomers");
   RunMonomers();
   MonoTimer.Stop();
   (*outfile)<<MonoTimer.PrintOut()<<std::endl;
   if(N>1){
      SmartTimer NMerTimer("Run NMers");
      RunNMers();
      NMerTimer.Stop();
      (*outfile)<<NMerTimer.PrintOut()<<std::endl;
   }
   Expansion<MBE> Expansion_(N);
   MBEProp<double> FinalEs=Expansion_.Property(*Frags_,*Energies_);
   (*outfile)<<FinalEs.PrintOut()<<std::endl;
   boost::shared_ptr<const LibParallel::Communicator> Comm=
         WorldComm->GetComm();
}

void LibFragDriver::RunCalc(int Start, int Stop) {
   typedef std::pair<int, SharedFrag> Pair_t;
   typedef MPITask<Pair_t> Task_t;
   std::vector<Task_t> Tasks;
   FragItr MonoI,MonoEnd;
   SmartTimer Overhead("Time to Make MPITasks");
   for (int i=Start; i<Stop; i++) {
      MonoEnd=Frags_->end(i);
      for (MonoI=Frags_->begin(i); MonoI!=MonoEnd; ++MonoI) {
         Pair_t temp(i, (*MonoI));
         Tasks.push_back(Task_t(temp, i+1));
      }
   }
   SmartTimer CommTimer("Communications");
   MPIJob<Pair_t> PMan(Tasks);
   PythonFxn<> RunCalc("wrappers", "new_run_calc", "s i");
   std::vector<double> TempEngy;
   Overhead.Stop();
   (*outfile)<<Overhead.PrintOut()<<std::endl;
   CommTimer.Start();
   for(Pair_t NMer=PMan.Begin();!PMan.Done();NMer=PMan.Next()) {
      CommTimer.Stop();
      std::string Frag=MakeMol(NMer.second);
      RunCalc(Frag.c_str(),1);
      double Egy=Process::environment.globals["CURRENT ENERGY"];
      TempEngy.push_back(Egy);
      CommTimer.Resume();
   }
   CommTimer.Stop();
   (*outfile)<<CommTimer.PrintOut()<<std::endl;
   SmartTimer SynchTimer("Sychronization time");
   std::vector<double> TempEgys2=PMan.Synch(TempEngy, 1);
   SynchTimer.Stop();
   (*outfile)<<SynchTimer.PrintOut()<<std::endl;
   for (int i=Start,counter=0; i<Stop; i++) {
      MonoEnd=Frags_->end(i);
      MonoI=Frags_->begin(i);
      for (int j=0;MonoI!=MonoEnd;++MonoI,++j)
         Energies_->Change(i,j,TempEgys2[counter++]);
   }
}

void LibFragDriver::RunNMers() {RunCalc(1, Frags_->N());}

std::string LibFragDriver::MakeMol(
      boost::shared_ptr<LibMolecule::Molecule> Mol) const {
   std::stringstream Molecule;
   Molecule<<Mol->PrintOut(0);
   Molecule<<"symmetry c1"<<std::endl;
   Molecule<<"no_reorient"<<std::endl;
   Molecule<<"no_com"<<std::endl;
   return Molecule.str();
}

}} //End namespaces
