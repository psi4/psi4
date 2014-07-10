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

#include "BasisSet.hpp"
#include "Molecule.hpp"
#include <sstream>
namespace JKFactory {

/**********************Primitive Member Functions**************************/
void Primitive::PrintOut()const{
	std::stringstream hp;
	hp<<GetBeta()<<" "<<GetC()<<std::endl;
	this->Print(hp.str());
}

/*********************Shell Member Functions*****************************/

void Shell::AddPrimitive(double c0,double beta0){
	Gi.push_back(pPrim(new Primitive(c0,beta0)));
}

std::string Shell::l2string(const int l)const{
	std::string symbol;
	switch(l){
		case (0): {symbol="S";break;}
		case (1): {symbol="P";break;}
		case (2): {symbol="D";break;}
		case (3): {symbol="F";break;}
		case (4): {symbol="G";break;}
		case (5): {symbol="H";break;}
		case (6): {symbol="I";break;}
		default: {
			Error("Angular momentum quantum number is too large");
			break;
		}
	}
	return symbol;
}

void Shell::PrintOut()const{
	std::stringstream hp;
	hp<<l2string(l)<<"   "<<GetNPrims()<<"   "<<1.00<<std::endl;
	for(int prim=0;prim<GetNPrims();prim++){
		(*this)[prim]->PrintOut();
	}
}

/**************BasisSet Member Functions******************************/
BasisSet::BasisSet()
		: NBasis(0), NShells(0), max_momentum(0), max_prims(0),
		  max_prims_ID(0), max_basis(0), NPrims(0){
}

void BasisSet::PrintOut()const{}

int BasisSet::GetNPrims(){
	if(NPrims==0) GetMaxPrim(max_prims_ID);
	return NPrims;
}
/**************AtomicBasisSet Member Functions***********************/
AtomicBasisSet::AtomicBasisSet(const double *Carts){
	for(int i=0;i<3;i++)carts[i]=Carts[i];
}

void AtomicBasisSet::PrintOut()const{
	for(int i=0;i<Shells.size();i++){
		(*this)[i]->PrintOut();
	}
}

void AtomicBasisSet::AddShell(int l,bool IsCart){
	if(l>max_momentum) max_momentum=l;
	if(IsCart) Shells.push_back(pShell(new CartShell(l)));
	else Shells.push_back(pShell(new SphereShell(l)));
	NShells++;
}
int AtomicBasisSet::GetMaxPrim(int &ID){
	if(max_prims==0){
		NPrims=0;
		for(int i=0;i<Shells.size();i++){
			int prims=(*this)[i]->GetNPrims();
			if(prims>max_prims){
				max_prims=prims;
				max_prims_ID=i;
			}
			NPrims+=prims;
		}
	}
	ID=max_prims_ID;
	return max_prims;
}

int AtomicBasisSet::GetMaxBasis(){
	if(max_basis==0){
		for(int i=0;i<GetNShells();i++){
			int nbasis=(*this)[i]->GetNBasis();
			if(nbasis>max_basis) max_basis=nbasis;
		}
	}
	return max_basis;
}

int AtomicBasisSet::GetNBasis(){
	if(NBasis==0){
		for(int i=0;i<Shells.size();i++)
			NBasis+=(*this)[i]->GetNBasis();
	}
	return NBasis;
}
/*********************AOBasisSet Member Functions**************************/

AOBasisSet::AOBasisSet(pMol molecule){
	//These both start at 0
	ShellsPerAtom.push_back(0);
	BFsPerShell.push_back(0);
	for(int i=0;i<molecule->NAtoms();i++){
		AOBasis.push_back(pABS(new AtomicBasisSet(&(*molecule)(i,0))));
	}
}

int AOBasisSet::AbsShellNum(int Atom,int Shell){
	if(ShellsPerAtom.size()==1){
		for(int i=0;i<this->GetNAtoms();i++){
			ShellsPerAtom.push_back(
					ShellsPerAtom[i]+(*this)[i]->GetNShells());
		}
	}
	return ShellsPerAtom[Atom]+Shell;
}

int AOBasisSet::AbsBasisNum(int Atom,int Shell,int Basis){
	int ShellN=AbsShellNum(Atom,Shell);
	if(BFsPerShell.size()==1){
		for(int i=0;i<this->GetNAtoms();i++){
			pABS ABS=(*this)[i];
			for(int j=0;j<ABS->GetNShells();j++){
				int currentindex=BFsPerShell.size()-1;
				BFsPerShell.push_back(
						BFsPerShell[currentindex]+(*ABS)[j]->GetNBasis());
			}
		}
	}
	return BFsPerShell[ShellN]+Basis;
}

void AOBasisSet::PrintOut()const{
	std::stringstream hp;
	hp<<"Basis is made of "<<NBasis<<" "<<
			((*(*this)[0])[0]->IsCart()?"Cart":"Pure")<<
			" basis functions"<<std::endl;
	for(int i=0;i<AOBasis.size();i++){
		hp<<i+1<<"    "<<0<<std::endl;
		(*this)[i]->PrintOut();
	}
}

int AOBasisSet::GetNShells(){
	if(NShells==0){
		for(int i=0;i<AOBasis.size();i++)
			NShells+=(*this)[i]->GetNShells();
	}
	return NShells;
}

int AOBasisSet::GetMaxL(){
	if(max_momentum==0){
		for(int i=0;i<AOBasis.size();i++){
			int l=(*this)[i]->GetMaxL();
			if(l>max_momentum) max_momentum=l;
		}
	}
	return max_momentum;
}

int AOBasisSet::GetMaxPrim(int &ID){
	if(max_prims==0){
		int tempNshells=0,tempID;
		NPrims=0;
		for(int i=0;i<AOBasis.size();i++){
			NPrims+=(*this)[i]->GetNPrims();
			int prims=(*this)[i]->GetMaxPrim(tempID);
			if(prims>max_prims){
				max_prims=prims;
				max_prims_ID=tempNshells+tempID;
			}
			tempNshells+=(*this)[i]->GetNShells();
		}
	}
	ID=max_prims_ID;
	return max_prims;
}

int AOBasisSet::GetMaxBasis(){
	if(max_basis==0){
		for(int i=0;i<GetNAtoms();i++){
			int nbasis=(*this)[i]->GetMaxBasis();
			if(nbasis>max_basis) max_basis=nbasis;
		}
	}
	return max_basis;
}

int AOBasisSet::GetNBasis(){
	if(NBasis==0){
		for(int i=0;i<AOBasis.size();i++){
			NBasis+=(*this)[i]->GetNBasis();
		}
	}
	return NBasis;
}

}//End namespace JKFactory

