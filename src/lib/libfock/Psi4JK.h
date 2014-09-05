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

#ifndef PSI4JK_HPP_
#define PSI4JK_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>
namespace psi{
class Matrix;
class BasisSet;
typedef boost::shared_ptr<psi::Matrix> SharedMatrix;
typedef boost::shared_ptr<psi::BasisSet> SharedPsiBasis;
}

#if HAVE_JK_FACTORY
#include "libJKFactory/src/Interface.hpp"
#include "libJKFactory/src/BasisReorderer.hpp"
namespace psi{


class Psi4JK:private JKFactory::Interface{
	private:
		void MakeMolecule(SharedPsiBasis &PsiBasis,pMol &JKMol);
		void MakeBasis(SharedPsiBasis &PsiBasis,pMyBasis &JKBasis);
		SharedMatrix Unitary;
		bool IsIdentity;
		std::vector<double> ScaleFacts;
	public:
		Psi4JK(SharedPsiBasis &PsiBasis);
		void UpdateDensity(SharedMatrix &Density,SharedMatrix& J,
              SharedMatrix& K);
		void UpdateDensity(std::vector<SharedMatrix>& Density,
		      std::vector<SharedMatrix>& J,
              std::vector<SharedMatrix>& K);
		bool Shared(){return false;}

};
class PsiReorderer:public JKFactory::BasisSetReorderer{
   protected:
      void YourCartConverter(int l, int k, int& a,int & b, int &c);
      void YourPureConverter(int l, int k, int& m);
};
}
#else
namespace psi{
class Psi4JK{
    public:
        Psi4JK(SharedPsiBasis &PsiBasis){}
        void UpdateDensity(SharedMatrix &Density,
              SharedMatrix &J,SharedMatrix &K){}
        void UpdateDensity(std::vector<SharedMatrix>& Density,
              std::vector<SharedMatrix>& J,
              std::vector<SharedMatrix>& K){}
        bool Shared(){return true;}
};
}
#endif

#endif /* PSI4JK_HPP_ */
