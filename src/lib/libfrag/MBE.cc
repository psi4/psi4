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

#include "MBE.h"
#include "MBEFrag.h"
#include "../../bin/psi4/psi4.h"
#include "../../../include/psi4-dec.h"
#include <iomanip>
namespace psi{
namespace LibFrag{
void MBE::MakeIntersections(std::vector<NMerSet>& Systems){
	int oldN=N;
	for(int i=2;i<oldN;i++){
		this->SetN(i);
		Systems.push_back(NMerSet());
		MakeNmers(Systems[0],Systems[Systems.size()-1]);
	}
	this->SetN(oldN);
}
double MBE::Energy(const std::vector<NMerSet>& Systems,
      const std::vector<boost::shared_ptr<double[]> >& Energies){
	double energy=0;
	//Total energy of each set of n-mers
	std::vector<double> En;
	std::vector<double> Egys;
	psi::outfile->Printf(
	      "\n************ MBE Energy Analysis ************\n\n");
	if(Systems.size()!=N&&N!=1)
	   throw psi::PSIEXCEPTION("The number of systems is not consistent with"
	         " the MBE truncation order....\n");
	for(int i=0;i<N;i++){
		En.push_back(0);//Initialize our vector
		int nfrags=Systems[i].size();
		for(int j=0;j<nfrags;j++){
			En[i]+=Energies[i][j];
		}
		energy=NBodyE(i+1,Systems[0].size(),&En[0]);
		Egys.push_back(energy);
		psi::outfile->Printf(
		      "%2d-body approximate energy: %16.15f (a.u.)\n",i+1,energy);
	}
	for(int i=1;i<N;i++){
	   double corr=(i==0?Egys[i]:Egys[i]-Egys[i-1]);
	   psi::outfile->Printf(
	         "%2d-body correction: %16.15f (a.u.)\n",i+1,corr);
	}
	return energy;
}

double MBE::NBodyE(const int N,const int nfrags,const double *DeltaEs){
	double energy=0;
	if(N!=nfrags){
		for(int i=0;i<N;i++)
		energy+=Coef(nfrags,N,i+1)*DeltaEs[i];
	}
	else energy=DeltaEs[N-1];
	return energy;
}

}}//End namespaces



