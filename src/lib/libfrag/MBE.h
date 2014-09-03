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

#ifndef MBE_H_
#define MBE_H_

#include "GMBE.h"
#include <boost/math/special_functions/binomial.hpp>
namespace psi{
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
		///Computes and returns the energy
		double Energy(const std::vector<MBEFragSet>& Systems,
		      const std::vector<boost::shared_ptr<double[]> >& Energies,
		      std::string& RealName);
		bool RunFrags()const{return true;}
		bool IsGMBE()const{return false;}
};
}}


#endif /* MBE_H_ */
