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

#include "MinimalInterface.h"
#include "libmints/basisset_parser.h"
#include "libmints/basisset.h"
#include "libmints/matrix.h"
#include "libmints/wavefunction.h"
#include "libmints/molecule.h"
#include "libmints/twobody.h"
#include "libmints/integral.h"
#include "psi4-dec.h"
extern "C" {
   #include "gtfock/libcint/CInt.h"
   #include "gtfock/pfock/pfock.h"
}


typedef boost::shared_ptr<psi::BasisSetParser> SharedParser;
typedef boost::shared_ptr<psi::BasisSet> SharedBasis;
typedef boost::shared_ptr<psi::Matrix> SharedMatrix;

//Helper fxns
void MakeBasis(BasisSet** GTBasis,SharedBasis PsiBasis);
void MyBlock(double** Block,SharedMatrix Matrix, int NPRow, int NPCol);
void BlockDims(int,int&,int&,int&,int&,int&);
void Gather(SharedMatrix,double*,int,int,int,int,int);

void psi::MinimalInterface::Vectorize(
      SharedMatrix Mat,
      void (psi::MinimalInterface::*fxn)(std::vector<SharedMatrix>&)){
   std::vector<SharedMatrix>temp(1);
   temp[0]=Mat;
   (this->*(fxn))(temp);
}

psi::MinimalInterface::~MinimalInterface(){
   PFock_destroy(PFock_);
   CInt_destroyBasisSet(GTBasis_);
}

psi::MinimalInterface::MinimalInterface(const int NMats,
      const bool AreSymm):NPRow_(1),NPCol_(1){
    psi::Options& options = psi::Process::environment.options;
    boost::shared_ptr<psi::BasisSet> primary =
    		psi::BasisSet::pyconstruct_orbital(psi::Process::environment.molecule(),
        "BASIS", options.get_str("BASIS"));
   MakeBasis(&GTBasis_,primary);
   double IntThresh=
         psi::Process::environment.options["INTS_TOLERANCE"].to_double();
   int NBlkFock=5;//No idea what this does, but I was told it equals 5
   PFock_create(GTBasis_,NPRow_,NPCol_,NBlkFock,IntThresh,
         NMats,AreSymm,&PFock_);
   std::cout<<"GTFock Made\n";

}

void psi::MinimalInterface::SetP(std::vector<SharedMatrix>& Ps){
   PFock_setNumDenMat(Ps.size(),PFock_);
   std::cout<<"Number of density matrices set to: "<<Ps.size()<<std::endl;
   int startr,endr,startc,endc,stride;
   BlockDims(Ps[0]->nrow(),startr,endr,startc,endc,stride);
   for(int i=0;i<Ps.size();i++){
      double* Buffer;
      Ps[0]->print_out();
      MyBlock(&Buffer,Ps[i],NPRow_,NPCol_);
      std::cout<<"("<<startr<<","<<startc<<") to  ("
    		   <<endr<<","<<endc<<")\n";
      PFock_putDenMat(startr,endr,startc,endc,stride,Buffer,i,PFock_);
   }
   PFock_commitDenMats(PFock_);
   std::cout<<"Density matrix committed\n";
   PFock_computeFock(GTBasis_,PFock_);
}

void psi::MinimalInterface::GenGetCall(
      std::vector<SharedMatrix>& JorK,
      int value){
   int startr,endr,startc,endc,stride;
   BlockDims(JorK[0]->nrow(),startr,endr,startc,endc,stride);
   double* Block=new double[(endr-startr)*(endc-startc)];
   for(int i=0;i<JorK.size();i++){
      PFock_getMat(PFock_,(PFockMatType_t)value,i,startr,
                  endr,startc,endc,stride,Block);
      Gather(JorK[i],Block,startr,endr,startc,endc,stride);
   }
   delete [] Block;
}


void psi::MinimalInterface::GetJ(std::vector<SharedMatrix>& Js){
   GenGetCall(Js,(int)PFOCK_MAT_TYPE_J);
}

void psi::MinimalInterface::GetK(std::vector<SharedMatrix>& Ks){
   GenGetCall(Ks,(int)PFOCK_MAT_TYPE_K);
}

void MyBlock(double **Buffer,SharedMatrix Matrix, int NPRow, int NPCol){
   Buffer[0]=&((*Matrix)(0,0));
}

void Gather(SharedMatrix Result,double *Block,int startr, int endr,
      int startc,int endc, int stride){
   int nrows=endr-startr+1;
   int ncols=endc-startc+1;
   if(!Result){
      Result=SharedMatrix(new psi::Matrix(nrows,ncols));
   }
   memcpy(&((*Result)(0,0)),Block,sizeof(double)*nrows*ncols);
}


void BlockDims(int NBasis,int& StartRow,int& EndRow,
               int& StartCol,int& EndCol,int& Stride){
   StartRow=0;
   StartCol=0;
   EndRow=NBasis-1;
   EndCol=NBasis-1;
   Stride=EndCol-StartCol+1;
}

void MakeBasis(BasisSet** GTBasis,SharedBasis PsiBasis){
   CInt_createBasisSet(GTBasis);
   boost::shared_ptr<psi::Molecule> PsiMol=PsiBasis->molecule();
   int NAtoms=PsiMol->natom();
   int NPrims=PsiBasis->nprimitive();
   int NShells=PsiBasis->nshell();
   int IsPure=(PsiBasis->has_puream()?0:1);
   //Carts,then atomic numbers
   std::vector<double> X(NAtoms),Y(NAtoms),Z(NAtoms),
                       Alpha(NPrims),CC(NPrims);
   std::vector<int> ShellsPerAtom(NAtoms),Zs(NAtoms),
         PrimsPerShell(NShells),L(NShells);
   for(int i=0;i<NAtoms;i++){
      Zs[i]=PsiMol->Z(i);
      X[i]=PsiMol->x(i);
      Y[i]=PsiMol->y(i);
      Z[i]=PsiMol->z(i);
      ShellsPerAtom[i]=PsiBasis->nshell_on_center(i);
   }
   for(int i=0,total=0;i<NShells;i++){
      PrimsPerShell[i]=PsiBasis->shell(i).nprimitive();
      L[i]=PsiBasis->shell(i).am();
      for(int j=0;j<PrimsPerShell[i];j++){
         Alpha[total]=PsiBasis->shell(i).exp(j);
         CC[total++]=PsiBasis->shell(i).coef(j);
      }
   }
   CInt_importBasisSet((*GTBasis), NAtoms, &Zs[0],
         &X[0], &Y[0], &Z[0], NPrims, NShells, IsPure,
         &ShellsPerAtom[0],&PrimsPerShell[0],&L[0],&CC[0],&Alpha[0]);
}


