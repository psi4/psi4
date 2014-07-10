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

#include "BasisReorderer.hpp"
#include "BasisSet.hpp"
namespace JKFactory {
void BasisSetReorderer::BuildT(pMyBasis& Basis) {
   int NBasis=Basis->GetNBasis();
   int NAtoms=Basis->GetNAtoms();
   T=new double[NBasis*NBasis];
   for (int i=0; i<NBasis*NBasis; i++)
      T[i]=0.0;
   for (int i=0,index=0; i<NAtoms; i++) {
      pABS ABS=(*Basis)[i];
      int NShells=ABS->GetNShells();
      for (int j=0; j<NShells; j++) {
         int L=(*ABS)[j]->GetL();
         int NOrbs=(*ABS)[j]->GetNBasis();
         for (int k=0; k<NOrbs; k++) {
            int JKValue=0;
            if ((*ABS)[j]->IsCart()) {
               int a=0,b=0,c=0;
               YourCartConverter(L, k, a, b, c);
               JKValue=JKFactoryCartConverter(L, a, b, c);
            }
            else {
               int m=0;
               YourPureConverter(L, k, m);
               JKValue=JKFactoryPureConverter(L,m);
            }
            if(k!=JKValue)IsIdentity=false;
            T[(index+k)*NBasis+(index+JKValue)]=1;
         }
         index+=NOrbs;
      }
   }
}

int BasisSetReorderer::JKFactoryCartConverter(int L, int a, int b, int c) {
   //Write out the first couple L's for speed
   int returnvalue=0;
   switch (L) {
      case (0): {
         returnvalue=0;
         break;
      }
      case (1): {
         if (a==1&&b==0&&c==0) returnvalue=0;
         else if (a==0&&b==1&&c==0) returnvalue=1;
         else if (a==0&&b==0&&c==1) returnvalue=2;
         break;
      }
      case (2): {
         if (a==2&&b==0&&c==0) returnvalue=0;
         else if (a==1&&b==1&&c==0) returnvalue=1;
         else if (a==1&&b==0&&c==1) returnvalue=2;
         else if (a==0&&b==2&&c==0) returnvalue=3;
         else if (a==0&&b==1&&c==1) returnvalue=4;
         else if (a==0&&b==0&&c==2) returnvalue=5;
         break;
      }
      case (3): {
         if (a==3&& b==0&&c==0) returnvalue=0;
         else if (a==2&&b==1&&c==0) returnvalue=1;
         else if (a==2&&b==0&&c==1) returnvalue=2;
         else if (a==1&&b==2&&c==0) returnvalue=3;
         else if (a==1&&b==1&&c==1) returnvalue=4;
         else if (a==1&&b==0&&c==2) returnvalue=5;
         else if (a==0&&b==3&&c==0) returnvalue=6;
         else if (a==0&&b==2&&c==1) returnvalue=7;
         else if (a==0&&b==1&&c==2) returnvalue=8;
         else if (a==0&&b==0&&c==3) returnvalue=9;
         break;
      }
      default:{
         int index=0;
         for(int i =0;i<=L&&index>=0;i++){
             for(int j=0;j<=i&&index>=0;j++){
                 int lx = L - i;
                 int ly = i - j;
                 int lz = j;
                 if(lx==a&&ly==b&&lz==c){
                    returnvalue=index;
                    index=-1;
                 }
                 else index++;
             }
         }
         break;
      }
   }
   return returnvalue;
}

int BasisSetReorderer::JKFactoryPureConverter(int L, int m){
   int returnvalue=0;
      switch (L) {
         case (0): {
            returnvalue=0;
            break;
         }
         case (1): {
            if (m==0) returnvalue=2;
            else if (m==1) returnvalue=0;
            else if (m==-1) returnvalue=1;
            break;
         }
         case (2):{
            if(m==0)returnvalue=2;
            else if(m==1)returnvalue=1;
            else if(m==-1)returnvalue=3;
            else if(m==2)returnvalue=0;
            else if(m==-2)returnvalue=4;
            break;
         }
         case(3):{
if(m==0)returnvalue= 0;
else if(m==1)returnvalue= 1;
else if(m==-1)returnvalue= 2;
else if(m==2)returnvalue= 3;
else if(m==-2)returnvalue= 5;
else if(m==3)returnvalue= 6;
else if(m==-3)returnvalue= 4;
            break;
         }
      }
   return returnvalue;
}
}

