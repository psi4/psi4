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
#ifndef SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_
#define SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_
#include <vector>
#include <boost/shared_ptr.hpp>

class BasisSet;
class PFock;
namespace psi{
class Matrix;
class TwoBodyAOInt;

namespace Interface{
extern std::vector<boost::shared_ptr<TwoBodyAOInt> > Ints;
}

class MinimalInterface{
   private:
      BasisSet* GTBasis_;
      PFock* PFock_;
      int NPRow_;
      int NPCol_;
      void GenGetCall(std::vector<boost::shared_ptr<Matrix> >& Mat,
            int Value);
      void Vectorize(boost::shared_ptr<Matrix> Mat,
         void (MinimalInterface::*fxn)(
               std::vector<boost::shared_ptr<Matrix> >&)
         );
   public:
      MinimalInterface(const int NMats=1,const bool AreSymm=true);
      ~MinimalInterface();
      void GetJ(std::vector<boost::shared_ptr<Matrix> >& Js);
      void GetJ(boost::shared_ptr<Matrix> J){
         Vectorize(J,&MinimalInterface::GetJ);
      }
      void GetK(std::vector<boost::shared_ptr<Matrix> >& Ks);
      void GetK(boost::shared_ptr<Matrix> K){
         Vectorize(K,&MinimalInterface::GetK);
      }
      void GetH(boost::shared_ptr<Matrix> H);
      void SetP(std::vector<boost::shared_ptr<Matrix> >& Ps);
      void SetP(boost::shared_ptr<Matrix> P){
         Vectorize(P,&MinimalInterface::SetP);
      }
};

}

#endif /* SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_ */
