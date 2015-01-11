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
#include "libmints/matrix.h"
#include "libmints/basisset.h"
#include "libmints/molecule.h"
#include "libmints/sieve.h"
#include "libparallel2/MPITask.h"
#include "libparallel2/MPIJob.h"
#include "libparallel2/Algorithms.h"
#include "libmints/integral.h"
#include "libmints/twobody.h"
#include "JKFactory.h"
#include "libqt/qt.h"
namespace psi {
typedef boost::shared_ptr<Matrix> SharedMatrix;
typedef boost::shared_ptr<BasisSetParser> SharedParser;
typedef boost::shared_ptr<BasisSet> SharedBasis;
typedef boost::shared_ptr<Molecule> SharedMol;
typedef boost::tuple<int,int,int,int> ShellQuartet;
typedef MPITask<ShellQuartet> Task;
typedef std::pair<int,int> IntPair;

//This fxn uses an ERI sieve to return a list of significant shell quartets
void GetQuarts(std::vector<Task>& Quarts, SharedBasis primary) {
   double IntThresh=
     psi::Process::environment.options["INTS_TOLERANCE"].to_double();
   boost::shared_ptr<ERISieve> sieve_ =
        boost::shared_ptr<ERISieve>(new ERISieve(primary, IntThresh));
   std::vector<IntPair> SigShellPairs;
   int NAtoms=primary->molecule()->natom();
   for (int Atom1=0; Atom1<NAtoms; Atom1++) {
      int ShellsAtom1=primary->nshell_on_center(Atom1);
      for (int Atom2=0; Atom2<NAtoms && Atom2<=Atom1; Atom2++) {
         int ShellsAtom2=primary->nshell_on_center(Atom2);
         for (int P=0; P<ShellsAtom1; P++) {
            int AbsP=primary->shell_on_center(Atom1, P);
            for (int Q=0; Q<ShellsAtom2; Q++) {
               int AbsQ=primary->shell_on_center(Atom2, Q);
               //This may need to be a continue if the shells aren't numbered
               //consecutively across consecutive atoms
               if (AbsQ>AbsP)break;
               if(sieve_->shell_pair_significant(AbsP, AbsQ))
                  SigShellPairs.push_back(IntPair(AbsP,AbsQ));
            }
         }
      }
   }

   //We will assume we found the significant shell pairs in lexical order
   std::vector<IntPair>::iterator Pair1,Pair2,done;
   done=SigShellPairs.end();
   for(Pair1=SigShellPairs.begin();Pair1!=done;++Pair1){
      int P=(*Pair1).first;
      int Q=(*Pair1).second;
      int Nbf=primary->shell(P).nfunction()*primary->shell(Q).nfunction();
      for(Pair2=SigShellPairs.begin();Pair2<=Pair1;++Pair2){
         int R=(*Pair2).first;
         int S=(*Pair2).second;
         int Nbf2=Nbf*primary->shell(R).nfunction()*primary->shell(S).nfunction();
         Quarts.push_back(Task(ShellQuartet(P,Q,R,S),Nbf2));
         //std::cout<<"Adding quartet: "<<P<<" "<<Q<<" "<<R<<" "<<S<<std::endl;
      }
   }

}

void JKFactory::BuildJandK(SharedMatrix Rho) {
   psi::Options& options=psi::Process::environment.options;
   SharedBasis primary=psi::BasisSet::pyconstruct_orbital(
         psi::Process::environment.molecule(), "BASIS",
         options.get_str("BASIS"));
   int NShell=primary->nshell();
   int NBasis=primary->nbf();
   int NBasis2=NBasis*NBasis;
   std::vector<int> AtomStarts,NBfPShell(NShell);

   //Get significant shell quartets as determined by ERI screening
   std::vector<Task> Quarts;
   timer_on("Make MPI Tasks");
   GetQuarts(Quarts,primary);
   timer_off("Make MPI Tasks");

   int MaxFxns=primary->max_function_per_shell();
   int NUniqueBuffers=6;
   std::vector<double>JBuffer(NBasis2,0);
   std::vector<double>KBuffer(NBasis2,0);

   boost::shared_ptr<IntegralFactory> IntFac(
         new IntegralFactory(primary,primary,primary,primary));
   boost::shared_ptr<TwoBodyAOInt> Integrals(IntFac->erd_eri());
   MPIJob<ShellQuartet> Job(Quarts);
   bool touched=false;
   timer_on("Main loop");
   for(ShellQuartet quart=Job.Begin();!Job.Done();quart=Job.Next()){
      int Shells[4];
      //Can't be in a loop because template param needs to be known at
      //compile time
      Shells[0]=quart.get<0>();Shells[1]=quart.get<1>();
      Shells[2]=quart.get<2>();Shells[3]=quart.get<3>();
      timer_on("Computing Integrals");
      int NInts=
            Integrals->compute_shell(Shells[0],Shells[1],Shells[2],Shells[3]);
      timer_off("Computing Integrals");
      if(NInts!=0){
         const double *IntBuffer=Integrals->buffer();
         int Sizes[4],Offsets[4];
         for(int i=0;i<4;i++){
            Sizes[i]=primary->shell(Shells[i]).nfunction();
            Offsets[i]=primary->shell_to_basis_function(Shells[i]);
         }
         double prefactor=1.0;
         if(Shells[0]==Shells[1])prefactor*=0.5;
         if(Shells[2]==Shells[3])prefactor*=0.5;
         if(Shells[0]==Shells[2]&&Shells[1]==Shells[3])prefactor*=0.5;
         timer_on("DGEMM");
         for(int p=0;p<Sizes[0];p++){
            int Absp=p+Offsets[0];
         for(int q=0;q<Sizes[1];q++){
            int Absq=q+Offsets[1];
         for(int r=0;r<Sizes[2];r++){
            int Absr=r+Offsets[2];
         for(int s=0;s<Sizes[3];s++){
            double IntVal=prefactor*(*IntBuffer);
            int Abss=s+Offsets[3];
            //J(p,q)+=D(r,s)*I(r,s,p,q)
            JBuffer[Absp*NBasis+Absq]+=
                  ((*Rho)(Absr,Abss)+(*Rho)(Abss,Absr))*IntVal;
            //J(r,s)+=D(p,q)*I(r,s,p,q)
            JBuffer[Absr*NBasis+Abss]+=
                  ((*Rho)(Absp,Absq)+(*Rho)(Absq,Absp))*IntVal;
            //K(p,r)+=D(q,s)*I(r,s,p,q)
            KBuffer[Absp*NBasis+Absr]+=(*Rho)(Absq,Abss)*IntVal;
            //K(p,s)+=D(q,r)*I(r,s,p,q)
            KBuffer[Absp*NBasis+Abss]+=(*Rho)(Absq,Absr)*IntVal;
            //K(q,r)+=D(p,s)*I(r,s,p,q)
            KBuffer[Absq*NBasis+Absr]+=(*Rho)(Absp,Abss)*IntVal;
            //K(q,s)+=D(p,r)*I(r,s,p,q)
            KBuffer[Absq*NBasis+Abss]+=(*Rho)(Absp,Absr)*IntVal;
            IntBuffer++;
         }}}}
         timer_off("DGEMM");
      }
   }
   timer_off("Main loop");
   timer_on("Cleanup");
   if(!J_)J_=SharedMatrix(new Matrix(NBasis,NBasis));
   if(!K_)K_=SharedMatrix(new Matrix(NBasis,NBasis));
   std::vector<double> MatBuffer;
   MatBuffer=Job.Reduce(JBuffer,NBasis2,LibParallel::ADD);
   ::memcpy(&(*J_)(0,0),&MatBuffer[0],NBasis2*sizeof(double));
   J_->scale(2.0);
   J_->hermitivitize();
   MatBuffer=Job.Reduce(KBuffer,NBasis2,LibParallel::ADD);
   ::memcpy(&(*K_)(0,0),&MatBuffer[0],NBasis2*sizeof(double));
   K_->scale(2.0);
   K_->hermitivitize();
   timer_off("Cleanup");
}

}

