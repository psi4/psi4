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

#include "Molecule.hpp"

namespace JKFactory{

Atom::Atom(const int Z_,const double* Carts_):
		Z(Z_),NElectrons(Z_){
	for(int comp=0;comp<3;comp++)Carts[comp]=Carts_[comp];
}

void Atom::PrintOut()const{
	std::stringstream message;
	message<<Z<<" ";
	for(int comp=0;comp<3;comp++)message<<Carts[comp]<<" ";
	message<<std::endl;
	this->Print(message.str());
}

Molecule::Molecule(const std::vector<int>& Zs,const std::vector<double>& Carts,
		const int charge):
	NElectrons(0){
	if((3*Zs.size())!=Carts.size())
		Error("Number of atomic numbers does not equal\
			1/3 the number of carts");
	for(int atom=0;atom<Zs.size();atom++){
		Atoms.push_back(pAtom(new Atom(Zs[atom],&Carts[3*atom])));
		NElectrons+=(*this)[atom]->NElec();
	}
	NElectrons-=charge;
}

void Molecule::PrintOut()const{
	for(int atom=0;atom<this->NAtoms();atom++){
		(*this)[atom]->PrintOut();
	}
}

//This is just a very long mapping from atomic numbers to
//strings, read at your own risk
std::string Atom::ZtoAt(const int Z) const{
	std::string At;
	switch(Z){
		case (0): {
			this->Print("You picked Ez, although a \
				valid choice in the Mass Effect series, it is\
				 not one in real life :-(" );
			break;
		}
		case (1): {At="H";break;}
		case (2): {At="HE";break;}
		case (3): {At="LI";break;}
		case (4): {At="BE";break;}
		case (5): {At="B";break;}
		case (6): {At="C";break;}
		case (7): {At="N";break;}
		case (8): {At="O";break;}
		case (9): {At="F";break;}
		case (10):{At="NE";break;}
		case (11): {At="NA";break;}
		case (12): {At="MG";break;}
		case (13): {At="AL";break;}
		case (14): {At="SI";break;}
		case (15): {At="P";	break;}
		case (16): {At="S";	break;}
		case (17): {At="CL";break;}
		case (18): {At="AR";break;}
		case (19): {At="K";break;}
		case (20): {At="CA";break;}
		case (21): {At="SC";break;}
		case (22): {At="TI";break;}
		case (23): {At="V";break;}
		case (24): {At="CR";break;}
		case (25): {At="MN";break;}
		case (26): {At="FE";break;}
		case (27): {At="CO";break;}
		case (28): {At="NI";break;}
		case (29): {At="CU";break;}
		case (30): {At="ZN";break;}
		case (31): {At="GA";break;}
		case (32): {At="GE";break;}
		case (33): {At="AS";break;}
		case (34): {At="SE";break;}
		case (35): {At="BE";break;}
		case (36): {At="KR";break;}
		default: {
			this->Print("I got tired of typing out symbols at Z=36,\
				so your symbol isn't coded yet. (Or you gave me a Z<0)\n");
		}
	}
	return At;
}

}

