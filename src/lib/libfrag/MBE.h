/*
 * MBE.h
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#ifndef MBE_H_
#define MBE_H_

#include "GMBE.h"
#include <boost/math/special_functions/binomial.hpp>
namespace LibFrag{
///Specialization of the GMBE to non-intersecting fragments
class MBE:public GMBE{
	private:
		///Returns N choose M
		double Binomial(const int N,const int M){return boost::math::binomial_coefficient<double>(N,M);}
		/** \brief Computes the total energy for an N-body expansion of nfrags
		 *
		 * \param[in] N Truncation order of the MBE
		 * \param[in] nfrags Number of fragments
		 * \param[in] Ens Array of summed energies, such that Ens[0]=sum of monomer energies, Ens[1]=sum of dimer energies, etc.
		 */
		double NBodyE(const int N,const int nfrags, const double *Ens);
		///Returns the coefficient for the m-th term of an N-body expansion of nfrags.  m=1 for Ens[0], m=N for Ens[N-1]
		double Coef(const int nfrags,const int N,const int m){return Phase(N-m)*Binomial(nfrags-m-1,N-m);}
	public:
		///Sets MBE truncation order to N, default is 1
		MBE(int N=1):GMBE(N){}
		///Function that computes the intersections of the NMers
		void MakeIntersections(std::vector<NMerSet>& Systems);
		///Computes and returns the energy
		double Energy(const std::vector<NMerSet>& Systems, const std::vector<double*>& Energies);
};
}


#endif /* MBE_H_ */
