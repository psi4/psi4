/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include "MinimalInterface.h"
#include "libmints/basisset_parser.h"
#include "libmints/basisset.h"
#include "libmints/matrix.h"
#include "psi4/src/lib/libmints/wavefunction.h"
#include "psi4/src/lib/libmints/molecule.h"
#include "libmints/twobody.h"
#include "libmints/integral.h"
#include "psi4/include/psi4-dec.h"
#include "../libparallel2/Communicator.h"
#include "../libparallel2/ParallelEnvironment.h"
#include "../libparallel2/Algorithms.h"
extern "C" {
   #include "CifiedFxns.h"
   #include "gtfock/libcint/CInt.h"
   #include "gtfock/pfock/pfock.h"
}


typedef boost::shared_ptr<psi::BasisSetParser> SharedParser;
typedef boost::shared_ptr<psi::BasisSet> SharedBasis;
typedef boost::shared_ptr<psi::Matrix> SharedMatrix;
typedef boost::shared_ptr<const psi::LibParallel::Communicator>
        ConstSharedComm;
typedef boost::shared_ptr<psi::LibParallel::Communicator>
        SharedComm;
//Helper fxns
void MakeBasis(BasisSet** GTBasis,SharedBasis PsiBasis);
void SplitProcs(int&,int&);

typedef void (psi::MinimalInterface::*GetPutFxn)(
                         std::vector<SharedMatrix>& );


void psi::MinimalInterface::Vectorize(
      SharedMatrix Mat,GetPutFxn fxn){
   std::vector<SharedMatrix>temp(1,Mat);
   (this->*(fxn))(temp);
}

psi::MinimalInterface::~MinimalInterface(){
   PFock_destroy(PFock_);
   CInt_destroyBasisSet(GTBasis_);
}

psi::MinimalInterface::MinimalInterface(const int NMats,
      const bool AreSymm):NPRow_(1),NPCol_(1),
                          StartRow_(0),StartCol_(0),
                          EndRow_(0),EndCol_(0),Stride_(0),NBasis_(0){
    SetUp();
    SplitProcs(NPRow_,NPCol_);
    psi::Options& options = psi::Process::environment.options;
    SharedBasis primary = psi::BasisSet::pyconstruct_orbital(
    		                  psi::Process::environment.legacy_molecule(),
                              "BASIS", options.get_str("BASIS"));
   NBasis_=primary->nbf();
   BlockDims(NBasis_);
   MakeBasis(&GTBasis_,primary);
   //I don't know how GTFock works
   double IntThresh=
         psi::Process::environment.options["INTS_TOLERANCE"].to_double();
   //It appears I can have GTFock figure this value out if I hand it
   //a negative value.
   int NBlkFock=-1;
   PFock_create(GTBasis_,NPRow_,NPCol_,NBlkFock,IntThresh,
         NMats,AreSymm,&PFock_);
}

void psi::MinimalInterface::SetP(std::vector<SharedMatrix>& Ps){
   PFock_setNumDenMat(Ps.size(),PFock_);
   double* Buffer;
   for(int i=0;i<Ps.size();i++){
      MyBlock(&Buffer,Ps[i]);
      PFock_putDenMat(StartRow_,EndRow_,
                      StartCol_,EndCol_,Stride_,Buffer,i,PFock_);
   }
   PFock_commitDenMats(PFock_);
   PFock_computeFock(GTBasis_,PFock_);
   delete [] Buffer;
}

void psi::MinimalInterface::GenGetCall(
      std::vector<SharedMatrix>& JorK,
      int value){
   int nrows=(EndRow_-StartRow_+1);
   int ncols=(EndCol_-StartCol_+1);
   double* Block=new double[nrows*ncols];
   for(int i=0;i<JorK.size();i++){
      memset(Block,0.0,sizeof(double)*nrows*ncols);
      PFock_getMat(PFock_,(PFockMatType_t)value,i,StartRow_,
                  EndRow_,StartCol_,EndCol_,Stride_,Block);
      Gather(JorK[i],Block);
   }
   delete [] Block;
}


void psi::MinimalInterface::GetJ(std::vector<SharedMatrix>& Js){
   GenGetCall(Js,(int)PFOCK_MAT_TYPE_J);
   for(int i=0;i<Js.size();i++)Js[i]->scale(0.5);
}

void psi::MinimalInterface::GetK(std::vector<SharedMatrix>& Ks){
   GenGetCall(Ks,(int)PFOCK_MAT_TYPE_K);
   for(int i=0;i<Ks.size();i++)Ks[i]->scale(-1.0);
}

void psi::MinimalInterface::MyBlock(double **Buffer,
                                    SharedMatrix Matrix){
   int nrows=EndRow_-StartRow_+1;
   (*Buffer)=new double[nrows*Stride_];
   for(int row=StartRow_;row<=EndRow_;++row){
      for(int col=StartCol_;col<=EndCol_;++col){
         (*Buffer)[(row-StartRow_)*Stride_+(col-StartCol_)]=
               (*Matrix)(row,col);
      }
   }
}

void FillMat(SharedMatrix Result, int RowStart, int RowEnd,
      int ColStart, int ColEnd, double* Buffer){
   int nrows=RowEnd-RowStart+1;
   int ncols=ColEnd-ColStart+1;
   for(int row=RowStart;row<=RowEnd;row++){
      for(int col=ColStart;col<=ColEnd;col++){
         (*Result)(row,col)=
               Buffer[(row-RowStart)*ncols+(col-ColStart)];
      }
   }
}


void psi::MinimalInterface::Gather(SharedMatrix Result,
                                   double *Block){
   if(!Result)Result=SharedMatrix(new psi::Matrix(NBasis_,NBasis_));
   SharedMatrix temp(new psi::Matrix(NBasis_,NBasis_));
   FillMat(temp,StartRow_,EndRow_,StartCol_,EndCol_,Block);
   ConstSharedComm Comm=psi::WorldComm->GetComm();
   Comm->AllReduce(&(*temp)(0,0),NBasis_*NBasis_,&(*Result)(0,0),LibParallel::ADD);
}

void psi::MinimalInterface::BlockDims(const int NBasis){
   ConstSharedComm Comm=psi::WorldComm->GetComm();
   int MyRank=Comm->Me();
   int IDs[2];
   //Divide the mat into NPRow_ by NPCol_ blocks
   //Note the following works because I know NPRow_*NPCol_=NProc
   IDs[0]=MyRank/NPCol_;//Row of my block
   IDs[1]=MyRank%NPCol_;//Col of my block
   for(int i=0;i<2;i++){
      SharedComm TempComm=Comm->MakeComm(IDs[i]);
      int MyDim=TempComm->Me();
      int & DimStart=(i==0?StartRow_:StartCol_);
      int & DimEnd=(i==0?EndRow_:EndCol_);
      int & NP=(i==0?NPRow_:NPCol_);
      int BFsPerDimPerBlock=NBasis/NP;
      DimStart=IDs[i]*BFsPerDimPerBlock;
      //This needs to be the last actually usable value
      DimEnd=(IDs[i]==(NP-1)?NBasis:DimStart+BFsPerDimPerBlock)-1;
      TempComm->FreeComm();
   }
   Stride_=EndCol_-StartCol_+1;
}

void SplitProcs(int& NPRow, int& NPCol){
   ConstSharedComm Comm=psi::WorldComm->GetComm();
   int NProc=Comm->NProc();
   NPRow=(int)floor(std::sqrt((double)NProc));
   bool done=false;
   while (!done) {
      if (NProc%NPRow==0) {
         NPCol=NProc/NPRow;
         done=true;
      }
      else NPRow--;
   }
}

void MakeBasis(BasisSet** GTBasis,SharedBasis PsiBasis){
   CInt_createBasisSet(GTBasis);
   boost::shared_ptr<psi::Molecule> PsiMol=PsiBasis->molecule();
   int NAtoms=PsiMol->natom();
   int NPrims=PsiBasis->nprimitive();
   int NShells=PsiBasis->nshell();
   int IsPure=(PsiBasis->has_puream()?1:0);
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


