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
#include <sstream>
#include <boost/python.hpp>
#include "psi4-dec.h"
#include "process.h"
#include "LibFragDriver.h"
#include "libmolecule/LibMolecule.h"
#include "../libmints/molecule.h"
#include "libparallel2/LibParallel2.h"
#include "libPsiUtil/PythonFxn.h"
#include "libPsiUtil/Timer.h"
#include "libPsiUtil/ProgressBar.h"
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
      if (SCellSides[0]+SCellSides[1]+SCellSides[2]>=1) {
         boost::shared_ptr<LibMolecule::SuperCell> temp(
               new LibMolecule::SuperCell((*UC),SCellSides));
         MySys=SharedFrags(
               new LibMolecule::FragmentedSystem(*temp, N));
      }
      else MyMol=UC;
   }
   if (!MySys)
      MySys=SharedFrags(new LibMolecule::FragmentedSystem(*MyMol, N));
   return MySys;
}

LibFragDriver::LibFragDriver(const std::string& MethodName){
   outfile->MakeBanner("Many-Body Expansion Module");
   Options options=psi::Process::environment.options;
   int N=options["TRUNCATION_ORDER"].to_integer();
   SharedMol AMol=psi::Process::environment.molecule();
   Frags_=MakeMolecule(N,AMol);
   RunMonomers(MethodName);
   if(N>1)
      RunNMers(MethodName);
   //Expansion<MBE> Expansion_(N);
   Expansion<CanonicalMBE> Expansion_(N);
   for(int i=0;i<Energies_.size();i++){
      MBEProp<double> FinalEs=Expansion_.Property(*Frags_,*(Energies_[i]));
      (*outfile)<<FinalEs.PrintOut()<<std::endl;
   }
}

typedef MPITask<SharedFrag> Task_t;
std::vector<Task_t> MakeTasks(const int Start, const int Stop,const SharedFrags& Frags_){
   std::vector<Task_t> Tasks;
   FragItr MonoI,MonoEnd;
   //SmartTimer Overhead("Time to Make MPITasks");
   for (int i=Start; i<Stop; i++) {
      MonoEnd=Frags_->end(i);
      for (MonoI=Frags_->begin(i); MonoI!=MonoEnd; ++MonoI)
         Tasks.push_back(Task_t((*MonoI), i+1));
   }
   //Overhead.Stop();
   //(*outfile)<<Overhead.PrintOut()<<std::endl;
   return Tasks;
}


void LibFragDriver::RunCalc(const std::string& MethodName,int Start, int Stop) {
   std::vector<Task_t> Tasks=MakeTasks(Start,Stop,Frags_);

   //Make a progress bar, by guessing each processor gets same num of tasks
   boost::shared_ptr<const LibParallel::Communicator> Comm=
         WorldComm->GetComm();
   long unsigned int Increment=floor(Tasks.size()/Comm->NProc());
   ProgressBar MyBar((Increment==0?1:Increment));

   PythonFxn<boost::python::dict> RunCalc("wrappers", "new_run_calc", "s s i");

   //Array for collecting energies, the number of energy types, and their names
   std::vector<double> TempEngy;
   int NEgys=0;
   boost::python::list* Keys=NULL;

   //Start Comm timer here since there are some comms in the MPIJob creation
   SmartTimer CommTimer("Communications"),WorkTimer("NMer Energy Computation");
   MPIJob<SharedFrag> PMan(Tasks);
   for(SharedFrag NMer=PMan.Begin();!PMan.Done();NMer=PMan.Next()) {
      CommTimer.Stop();
      boost::python::dict Egys=
            RunCalc(MethodName.c_str(),MakeMol(NMer).c_str(),1);
      ++MyBar;
      if(NEgys==0){
         Keys=new boost::python::list(Egys.keys());
         NEgys=boost::python::len((*Keys));
      }
      for(int i=0;i<NEgys;i++)
         TempEngy.push_back(boost::python::extract<double>(Egys[(*Keys)[i]]));
      //TempEngy.push_back(Process::environment.globals["CURRENT ENERGY"]);
      CommTimer.Resume();
   }
   CommTimer.Stop();
   WorkTimer.Stop();
   (*outfile)<<std::endl<<CommTimer.PrintOut()<<std::endl;
   (*outfile)<<WorkTimer.PrintOut()<<std::endl;
   if(Energies_.size()==0){
      for(int i=0;i<NEgys;i++){
         const std::string EgyName=
               boost::python::extract<std::string>((*Keys)[i]);
         Energies_.push_back(
               boost::shared_ptr<MBEProp<double> >(
                     new MBEProp<double>(Frags_->N(),EgyName)));
      }
      //Energies_.push_back(
      //      boost::shared_ptr<MBEProp<double> >(
      //            new MBEProp<double>(Frags_->N(),"CURRENT ENERGY")));
   }
   delete Keys;
   //SmartTimer SynchTimer("Sychronization time");
   std::vector<double> TempEgys2=PMan.Synch(TempEngy, NEgys);
   //SynchTimer.Stop();
   //(*outfile)<<SynchTimer.PrintOut()<<std::endl;

   std::vector<Task_t>::iterator TaskI=Tasks.begin(),TaskIEnd=Tasks.end();
   for(int counter=0;TaskI!=TaskIEnd;++TaskI){
      const LibMolecule::SerialNumber& SN=TaskI->GetLabel()->GetSN();
      for(int i=0;i<NEgys;i++)
         Energies_[i]->Change(SN.size()-1,SN,TempEgys2[counter++]);
   }
}

void LibFragDriver::RunNMers(const std::string& MethodName) {
   RunCalc(MethodName,1, Frags_->N());
}

void LibFragDriver::RunMonomers(const std::string& MethodName){
   RunCalc(MethodName,0,1);
}

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
