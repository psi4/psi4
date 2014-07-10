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
#ifndef BASISREORDERER_HPP_
#define BASISREORDERER_HPP_

#include "APIs/Interface.hpp"
namespace JKFactory{
class BasisSetReorderer{
         protected:
            double* T;
            int JKFactoryCartConverter(int l,int a, int b, int c);
            int JKFactoryPureConverter(int l,int m);
            virtual void YourCartConverter(int l, int k,int& a,int& b,
                  int& c)=0;
            virtual void YourPureConverter(int l, int k,int& m)=0;

         public:
            BasisSetReorderer():T(NULL),IsIdentity(true){}
            operator double*()const{return T;}
            void BuildT(pMyBasis& Basis);
            virtual ~BasisSetReorderer(){delete [] T;}
            bool IsIdentity;

   };

}


#endif /* BASISREORDERER_HPP_ */
