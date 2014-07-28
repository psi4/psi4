/*
 * GMBE.h
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#ifndef GMBE_H_
#define GMBE_H_
#include "LibFragTypes.h"
namespace LibFrag{

class GMBE{
   private:
       ///The multiplicity of each n-mer and intersection in the positive ints
       std::vector<int> PMults;
       ///The multiplicity of each intersection in the negative ints
       std::vector<int> NMults;
   protected:
		///The order of the GMBE (e.g. N=1, is a one-body expansion etc.)
		int N;
		///Returns -1^i, the phase of many terms in the GMBE
		double Phase(const int i){return (i%2==0?1:-1);}
	public:
		///Makes an N-body GMBE, N defaults to 1
		GMBE(int newN=1):N(newN){}
		///Allows you to change N
		void SetN(const int newN){N=newN;}
		///Given a set of Monomers, returns all possible unions of N monomers
		void MakeNmers(const NMerSet& Monomers,NMerSet& NMers);
		///Given a set of Monomers and Nmers, makes the intersections of the Nmers
		virtual void MakeIntersections(std::vector<NMerSet>& Systems);
		///Returns the total energy of this GMBE expansion
		virtual double Energy(const std::vector<NMerSet>& Systems,const std::vector<double*>& Energies);
		///No memory allocated
		virtual ~GMBE(){}
		/** \brief Returns true if we need to run the fragments
		 *
		 *  For the GMBE there is currently no known form that allows us
		 *  to recover the one-body energy for a n=2 or higher expansion.
		 *  Consequentially there is no reason to run the fragments, unless
		 *  we are performing electrostatic embedding.  When the MBE derives
		 *  from this class it will overwrite this behavior because it can
		 *  recover the one-body energy for arbitrary truncations.
		 *
		 */
		virtual bool RunFrags()const{return (N==1);}
		virtual bool IsGMBE()const{return true;}
};

}//End namespace



#endif /* GMBE_H_ */
