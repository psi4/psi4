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
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "Interface.hpp"
#include "Matrix.hpp"
#include "BasisSet.hpp"
#include "Molecule.hpp"
namespace JKFactory{

//This fxn makes Vec the size of NewMats and puts NewMats[0] in Vec[0], etc.
void UpdateVector(std::vector<double*>&  NewMats,
		std::vector<pMatrix>& Vec,const int dim){
        int CurrSize=Vec.size(),NewSize=NewMats.size();
		if(CurrSize<NewSize){
			int diff=NewSize-CurrSize;//Number of matrices they differ by
			for(int i=0;i<diff;i++){
			    //Need all matrices from CurrSize on
			    double* Mat=NewMats[CurrSize+i];
				Vec.push_back(pMatrix(new JKFactory::Matrix(Mat,dim)));
			}
		}
		//Update matrices up to CurrSize
		for(int i=0;i<CurrSize;i++)
			Vec[i]->Update(NewMats[i]);
}

void Interface::UpdateDensity(double* NewDensity){
	std::vector<double*> Temp;
	Temp.push_back(NewDensity);
	Interface::UpdateDensity(Temp);
}

void Interface::UpdateDensity(std::vector<double*>& NewDensitys){
	UpdateVector(NewDensitys,Densities,Basis->GetNBasis());
	BuildJK(Basis,System);
}

void Interface::UpdateJ(double* NewJ){
    std::vector<double*> Temp;
	Temp.push_back(NewJ);
	Interface::UpdateJ(Temp);
}

void Interface::UpdateJ(std::vector<double*>& NewJs){
   int nbasis=Basis->GetNBasis();
   for(int j=0;j<NewJs.size();j++){
      for(int i=0;i<nbasis*nbasis;i++)NewJs[j][i]=0.0;
   }
   UpdateVector(NewJs,Js,nbasis);
}

void Interface::UpdateK(double* NewK){
	std::vector<double*> Temp;
	Temp.push_back(NewK);
	Interface::UpdateK(Temp);
}

void Interface::UpdateK(std::vector<double*>& NewKs){
   int nbasis=Basis->GetNBasis();
   for(int j=0;j<NewKs.size();j++){
      for(int i=0;i<nbasis*nbasis;i++)NewKs[j][i]=0.0;
   }
   UpdateVector(NewKs,Ks,nbasis);
}

void Interface::PrintOut()const{
	std::stringstream hp;
	hp<<"Integral threshold: "<<IntThresh<<" a.u."<<std::endl;
	hp<<"All matrices are assumed to be ";
	if(!AreSymm)hp<<"NOT ";
	hp<<"symmetric."<<std::endl;
	this->Print(hp.str());
	System->PrintOut();
	Basis->PrintOut();
	GTFockInter::PrintOut();
}





}


