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
#ifndef GTFOCKINTER_HPP_
#define GTFOCKINTER_HPP_

#include "JKFactory.hpp"
#include <vector>
#include <boost/shared_ptr.hpp>
namespace JKFactory{
	class AOBasisSet;
	class Molecule;
	class Matrix;
}
typedef boost::shared_ptr<JKFactory::Molecule> pMol;
typedef boost::shared_ptr<JKFactory::AOBasisSet> pMyBasis;
typedef boost::shared_ptr<JKFactory::Matrix> pMatrix;
class BasisSet;
class PFock;
class GTFockInter:public JKFactory::JKFactoryBase{
	protected:
		///These are the densities that we are using to build J and K
		std::vector<pMatrix> Densities;
		///The resulting J matrices, note Js[0] came from Densities[0],etc.
		std::vector<pMatrix> Js;
		///The resulting K matrices, note Ks[0] came from Densities[0],etc.
		std::vector<pMatrix> Ks;
		///Actual call to tell GTFock to build J and K
		void BuildJK(const pMyBasis& BasisSet,const pMol& System);
		///The integral screening threshold
		double IntThresh;
		///Whether the density matrices are symmetric or not
		bool AreSymm;
	private:
		//No shared pointers, cause they are "C" classes
		BasisSet* GTBasis;
		PFock* GTFock;
		///Copies our basis set object into the GTFock one
		void CopyBasis(BasisSet* GTBasis,const pMyBasis& BasisSet,
				const pMol& Molecule);
		///Boots this interface up
		void Initialize(const pMyBasis& BasisSet, const pMol& Molecule);
	public:
		GTFockInter(double ScreenThresh,bool Symm=true):
			IntThresh(ScreenThresh),AreSymm(Symm),GTBasis(NULL),GTFock(NULL){}
		///Frees up memory associated with GTFock
		~GTFockInter();
		virtual void PrintOut()const;
};



#endif /* GTFOCKINTER_HPP_ */
