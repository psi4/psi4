/*
 * GMBE.h
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#ifndef GMBE_H_
#define GMBE_H_
#include <vector>
#include <boost/shared_ptr.hpp>
namespace LibFrag{
class MBEFrag;
class GMBE{
	protected:
		typedef boost::shared_ptr<MBEFrag> SharedFrag;
		typedef std::vector<SharedFrag> NMerSet;
		///The order of the GMBE (e.g. N=1, is a one-body expansion etc.)
		int N;
		///Returns -1^i, the phase of many terms in the GMBE
		double Phase(const int i){return (i%2==0?1:-1);}
	public:
		///Makes an N-body GMBE, N defaults to 1
		GMBE(int N_=1):N(N_){}
		///Allows you to change N
		void SetN(const int newN){N=newN;}
		///Given a set of Monomers, returns all possible unions of N monomers
		void MakeNmers(const NMerSet& Monomers,NMerSet& NMers);
		///Given a set of Monomers and Nmers, makes the intersections of the Nmers
		virtual void MakeIntersections(std::vector<NMerSet>& Systems)=0;
		///Returns the total energy of this GMBE expansion
		virtual double Energy(const std::vector<NMerSet>& Systems, const std::vector<double*>& Energies)=0;
		///No memory allocated
		virtual ~GMBE(){}
};

}



#endif /* GMBE_H_ */
