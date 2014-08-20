/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
extern "C"{
	#include "CInt.h"
    #include "cint_type.h"
	#include "PFock.h"
	#include "pfock_type.h"
}
#include "GTFockInter.hpp"
#include "BasisSet.hpp"
#include "Matrix.hpp"
#include "Molecule.hpp"
#include "MPIManager.hpp"
void GTFockInter::CopyBasis(BasisSet* GTbasis,const pMyBasis& basis,
		const pMol& molecule){
	GTbasis->natoms=basis->GetNAtoms();
	GTbasis->eid=(int *)malloc(sizeof(int)*GTbasis->natoms);
	GTbasis->xn=(double *)malloc(sizeof(double)*GTbasis->natoms);
	GTbasis->yn=(double *)malloc(sizeof(double)*GTbasis->natoms);
	GTbasis->zn=(double *)malloc(sizeof(double)*GTbasis->natoms);
	GTbasis->ncharge=(double *)malloc(sizeof(double)*GTbasis->natoms);
	GTbasis->nelectrons=molecule->NElec();
	GTbasis->basistype=((*(*basis)[0])[0]->IsCart()?0:1);
	GTbasis->lenatom0=basis->GetNAtoms();
	GTbasis->eptr=(int *)malloc(sizeof(int)*(GTbasis->natoms+1));
	GTbasis->lenshell0=basis->GetNShells();
	GTbasis->atom_start=(int *)malloc(sizeof(int)*(GTbasis->natoms+1));
	GTbasis->atom_start[GTbasis->natoms]=GTbasis->lenshell0;
	GTbasis->totnexp=basis->GetNPrims();
	GTbasis->ptrshell=(int*)malloc(sizeof(int)*(GTbasis->lenshell0+1));
	GTbasis->ptrshell[GTbasis->lenshell0]=GTbasis->totnexp;
	GTbasis->nexp0=(int*)malloc(sizeof(int)*GTbasis->totnexp);
	GTbasis->cc0=(double *)malloc(sizeof(double)*(GTbasis->totnexp));
	GTbasis->exp0=(double *)malloc(sizeof(double)*(GTbasis->totnexp));
	GTbasis->momentum0=(int *)malloc(sizeof(int)*GTbasis->lenshell0);
	GTbasis->nshells=basis->GetNShells();
	GTbasis->nfunctions=basis->GetNBasis();
	GTbasis->f_start_id=(int *)malloc(sizeof(int)*GTbasis->nshells);
	GTbasis->f_end_id=(int *)malloc(sizeof(int)*GTbasis->nshells);
	GTbasis->s_start_id=(int *)malloc(sizeof(int)*(GTbasis->natoms+1));
	GTbasis->s_start_id[GTbasis->natoms]=GTbasis->nshells;
	GTbasis->nexp=(int *)malloc(sizeof(int)*GTbasis->nshells);
	GTbasis->exp=(double **)malloc(sizeof(double *)*GTbasis->nshells);
	GTbasis->cc=(double **)malloc(sizeof(double *)*GTbasis->nshells);
	GTbasis->momentum=(int *)malloc(sizeof(int)*GTbasis->nshells);
	GTbasis->x=(double *)malloc(sizeof(double)*GTbasis->nshells);
	GTbasis->y=(double *)malloc(sizeof(double)*GTbasis->nshells);
	GTbasis->z=(double *)malloc(sizeof(double)*GTbasis->nshells);
	GTbasis->maxdim=basis->GetMaxBasis();
	GTbasis->max_momentum=basis->GetMaxL();
	GTbasis->max_nexp=basis->GetMaxPrim(GTbasis->max_nexp_id);
	int shelln=0,basisn=0,primn=0;
	for(int i=0;i<basis->GetNAtoms();i++){
		int Z=(*molecule)[i]->GetZ();
		GTbasis->eid[i]=Z;
		GTbasis->xn[i]=(*molecule)(i,0);
		GTbasis->yn[i]=(*molecule)(i,1);
		GTbasis->zn[i]=(*molecule)(i,2);
		GTbasis->ncharge[i]=(double)Z;
		GTbasis->eptr[Z]=i;
		GTbasis->atom_start[i]=shelln;
		boost::shared_ptr<JKFactory::AtomicBasisSet> ABS=(*basis)[i];
		GTbasis->s_start_id[i]=shelln;
		for(int j=0;j<ABS->GetNShells();j++){
			GTbasis->f_start_id[shelln]=basisn;
			boost::shared_ptr<JKFactory::Shell> shell=(*ABS)[j];
			basisn+=shell->GetNBasis();
			GTbasis->f_end_id[shelln]=basisn-1;
			GTbasis->x[shelln]=(*molecule)(i,0);
			GTbasis->y[shelln]=(*molecule)(i,1);
			GTbasis->z[shelln]=(*molecule)(i,2);
			GTbasis->nexp0[shelln]=shell->GetNPrims();
			GTbasis->nexp[shelln]=shell->GetNPrims();
			GTbasis->momentum[shelln]=shell->GetL();
			GTbasis->momentum[shelln]=shell->GetL();
			GTbasis->ptrshell[shelln]=primn;
			GTbasis->cc[shelln]=&(GTbasis->cc0[primn]);
			GTbasis->exp[shelln]=&(GTbasis->exp0[primn]);
			for(int k=0;k<shell->GetNPrims();k++){
				boost::shared_ptr<JKFactory::Primitive> prim=(*shell)[k];
				GTbasis->cc0[primn]=prim->GetC();
				GTbasis->exp0[primn]=prim->GetBeta();
				primn++;
			}
			shelln++;
		}
	}
}

void GTFockInter::Initialize(const pMyBasis& BasisSet,const pMol& System){
	std::string Comm=MPI->Comm();
	CInt_createBasisSet(&GTBasis);
	CopyBasis(GTBasis,BasisSet,System);
	MPI->Sync(Comm);
	int nprow=Js[0]->MPIData()->NPRow;
	int npcol=Js[0]->MPIData()->NPCol;
	int NMats=Js.size();
	int heapsize=0,stacksize=0;
	PFock_GAInit(GTBasis->nfunctions,nprow,npcol,NMats,heapsize,stacksize);
	MPI->Sync(Comm);
	int NBlkFock=5;
	PFock_create(GTBasis,nprow,npcol,NBlkFock,NMats,AreSymm,IntThresh,&GTFock);
}

void GTFockInter::BuildJK(const pMyBasis& BasisSet,const pMol& System){
	if(GTBasis==NULL && GTFock==NULL)Initialize(BasisSet,System);
	int startr=Js[0]->MPIData()->startrow;
	int endr=Js[0]->MPIData()->endrow;
	int startc=Js[0]->MPIData()->startcol;
	int endc=Js[0]->MPIData()->endcol;
	int stride=endc-startc+1;
	GTFock->committed=0;
	PFock_setNumDenMat(GTFock,Js.size());
	for(int i=0;i<Js.size();i++){
		double *MyBlock=Densities[i]->GetMyBlock();
		PFock_putDenMat(GTFock,i,startr,endr,startc,endc,MyBlock,stride);
	}
	PFock_commitDenMats(GTFock);
	PFock_computeFock(GTFock,GTBasis);
	for(int i=0;i<Js.size();i++){
		double* JBlock=&((Js[i]->MPIData())->MyBuffer[0]);
		PFock_getMat(GTFock,i,PFOCK_MAT_TYPE_J,startr,endr,startc,endc,
				JBlock,stride);

		Js[i]->Gather();
		if(Ks.size()>0){
			double* KBlock=&((Ks[i]->MPIData())->MyBuffer[0]);
			PFock_getMat(GTFock,i,PFOCK_MAT_TYPE_K,startr,endr,startc,endc,
				KBlock,stride);
			Ks[i]->Gather();
		}
	}
}


void GTFockInter::PrintOut()const{
	for(int i=0;i<Densities.size();i++){
		Densities[i]->PrintOut();
		Js[i]->PrintOut();
		Ks[i]->PrintOut();
	}
}

GTFockInter::~GTFockInter(){
	PFock_destroy(GTFock);
	CInt_destroyBasisSet(GTBasis);
}

