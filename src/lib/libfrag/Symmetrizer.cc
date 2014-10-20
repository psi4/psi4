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

#include "Symmetrizer.h"
#include "MBEFrag.h"
#include "libmints/vector.h"
#include "libmints/matrix.h"

namespace psi{
namespace LibFrag{

typedef boost::shared_ptr<const MBEFrag> cSharedFrag;



typedef boost::shared_ptr<Vector> SharedVec;
typedef std::vector<SharedVec> stdVecSharedVec;

inline void Initialize(boost::shared_ptr<stdVecSharedVec>& Carts,
      const CartSet<SharedAtom>& Frag){
   boost::shared_ptr<stdVecSharedVec> temp(new stdVecSharedVec);
   Carts=temp;
   int natoms=Frag.size();
   Matrix Itensor(*Frag.MoI());
   Matrix MoI(Itensor);
   Vector3 CoM(Frag.CoM());
   Vector EigenValues(3);
   Itensor.diagonalize(&MoI,&EigenValues);
   double tol=1.0e-4;
   //Eigenvalues come back in order, so 1!=3, unless 1==2 && 2==3
   bool e12=(abs(EigenValues[0]-EigenValues[1])<tol);
   bool e23=(abs(EigenValues[1]-EigenValues[2])<tol);
   if(e12||e23){
      throw PSIEXCEPTION("You actually don't have an asymmetric top...I guess I need to do more coding\n");
   }
   for(int i=0;i<natoms;i++){
      Carts->push_back(SharedVector(new Vector(3)));
      Vector OldCarts(3);
      for(int j=0;j<3;j++)OldCarts[j]=Frag.Object(Frag[i])->Carts()[j]-CoM[j];
      (*Carts)[i]->gemv(true,1.0,&MoI,&OldCarts,0.0);
   }
}


void Symmetrizer::RemoveDuplicates(MBEFragSet& FragSet)const{
   //Assume all Frags are same order
   int N=FragSet[0]->GetMBEOrder();
   int Size=FragSet.size();
   boost::shared_ptr<bool[]> Unique(new bool[Size]);
   std::vector<boost::shared_ptr<stdVecSharedVec> >Carts(Size);
   for(int fragi=0;fragi<Size;fragi++)Unique[fragi]=true;
   for(int fragi=0;fragi<Size;fragi++){
      if(Unique[fragi]){
         const CartSet<SharedAtom>& Fragi=FragSet[fragi]->Atoms();
         if(!Carts[fragi])Initialize(Carts[fragi],Fragi);
         int iSize=Fragi.size();
         for(int fragj=fragi+1;fragj<Size;fragj++){
            if(Unique[fragj]){
               const CartSet<SharedAtom>& Fragj=FragSet[fragj]->Atoms();
               if(!Carts[fragj])Initialize(Carts[fragj],Fragj);
               int jSize=Fragj.size();
               bool AreSame=(iSize==jSize);
               if(AreSame){
                  for(int atom=0;atom<iSize&&AreSame;atom++){
                     AreSame=false;
                     int Z1=Fragi.Object(Fragi[atom])->Z();
                     double q[3]={(*(*Carts[fragi])[atom])[0],(*(*Carts[fragi])[atom])[1],(*(*Carts[fragi])[atom])[2]};
                     for(int atom2=0;atom2<iSize&&!AreSame;atom2++){//Atom order isn't necessarily the same
                        int Z2=Fragj.Object(Fragj[atom2])->Z();
                        if(Z1==Z2){
                           int coordssame=0;
                           for(int i=0;i<3;i++){
                              double q2=(*(*Carts[fragj])[atom2])[i];
                              double diff=fabs(q[i]-q2);
                              if(diff<1e-1)coordssame++;
                           }
                           if(coordssame==3)AreSame=true;
                        }
                     }
                  }
               }
               if(AreSame){
                  Unique[fragj]=false;
                  ++(*FragSet.Frags_[fragi]);
               }
            }
         }
      }
   }

   std::vector<boost::shared_ptr<MBEFrag> > TempFrags;
   for(int i=0;i<Size;i++){
      if(Unique[i])TempFrags.push_back(FragSet.Frags_[i]);
   }
   FragSet.Frags_=TempFrags;
}



}}

