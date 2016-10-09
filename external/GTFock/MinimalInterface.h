/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_
#define SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_
#include <vector>
 #include "psi4/include/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <boost/shared_ptr.hpp>
 PRAGMA_WARNING_POP

class BasisSet;
class PFock;
namespace psi{
class Matrix;
class TwoBodyAOInt;


#ifdef HAVE_JK_FACTORY
namespace Interface{
extern std::vector<boost::shared_ptr<TwoBodyAOInt> > Ints;
}

class MinimalInterface{
   private:
      BasisSet* GTBasis_;
      PFock* PFock_;
      int NPRow_;
      int NPCol_;
      int StartRow_;
      int StartCol_;
      int EndRow_;
      int EndCol_;
      int Stride_;
      int NBasis_;
      void BlockDims(const int NBasis);
      void MyBlock(double **Buffer,boost::shared_ptr<Matrix> Matrix);
      void Gather(boost::shared_ptr<Matrix>,double*);
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
      void SetP(std::vector<boost::shared_ptr<Matrix> >& Ps);
      void SetP(boost::shared_ptr<Matrix> P){
         Vectorize(P,&MinimalInterface::SetP);
      }
      bool UseJKFactory()const{return true;}
};
#else
class MinimalInterface{
   public:
      MinimalInterface(const int NMats=1,const bool AreSymm=true){}
      ~MinimalInterface(){}
      void GetJ(std::vector<boost::shared_ptr<Matrix> >& Js){}
      void GetJ(boost::shared_ptr<Matrix> J){}
      void GetK(std::vector<boost::shared_ptr<Matrix> >& Ks){}
      void GetK(boost::shared_ptr<Matrix> K){}
      void SetP(std::vector<boost::shared_ptr<Matrix> >& Ps){}
      void SetP(boost::shared_ptr<Matrix> P){}
      bool UseJKFactory()const{return false;}
};
#endif
}//End psi namespace

#endif /* SRC_LIB_LIBJKFACTORY_SRC_MINIMALINTERFACE_H_ */
