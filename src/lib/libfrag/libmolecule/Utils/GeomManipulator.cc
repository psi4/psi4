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


#include "../Utils/GeomManipulator.h"
#include "MolItr.h"
#include "LibFragMolecule.h"
namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<double[]> SharedDouble;

double GeomManipulator::operator()(const int i, const int j)const{
   return Carts_[i*3+j];
}

SharedDouble GeomManipulator::GetCarts()const{
   SharedDouble DaCarts(new double[3*Mol_->NAtoms()]);
   MolItr CurrAtom=Mol_->Begin(),LastAtom=Mol_->End();
   for(int index=0;CurrAtom!=LastAtom;++CurrAtom)
      for(int i=0;i<3;i++)DaCarts[index++]=(*(*CurrAtom))[i];
   return DaCarts;
}

void GeomManipulator::SetCarts(boost::shared_ptr<double[]> NewCarts){
   MolItr CurrAtom=Mol_->Begin(),LastAtom=Mol_->End();
   for(int index=0;CurrAtom!=LastAtom;++CurrAtom){
      boost::shared_ptr<Atom> MyAtom=
            boost::const_pointer_cast<Atom>(*CurrAtom);
      for(int i=0;i<3;i++)(*MyAtom)[i]=NewCarts[index++];
   }
}

GeomManipulator::GeomManipulator(Molecule* Mol):
      Mol_(Mol){
   Carts_=GetCarts();
}

void GeomManipulator::Rotate(const double* RotMat,const bool trans){
   int NAtoms=Mol_->NAtoms();
   boost::shared_ptr<double[]> Result(new double [3*NAtoms]);
   DGEMM_Int(Result.get(),Carts_.get(),RotMat,NAtoms,3,3,false,trans);
   Carts_=Result;
}

void GeomManipulator::Scale(const double ScaleFac){
   for(int i=0;i<3*Mol_->NAtoms();i++)Carts_[i]*=ScaleFac;
}

void GeomManipulator::Translate(const double* TransVec,
      const double ScaleFac){
   int NAtoms=Mol_->NAtoms();
   double TempVec[3];
   for(int i=0;i<NAtoms;i++)
      Translate(TransVec,i,ScaleFac);
}
///Rotations around the Z,Y,X planes

void RotMat(double Rot[],const double Theta,
            const Plane RotPlane,const bool ClockWise=true){
   double CosA=cos(Theta),SinA=sin(Theta);
   Rot[0]=(RotPlane==XY||RotPlane==XZ?CosA:1);
   Rot[4]=(RotPlane==YZ||RotPlane==XY?CosA:1);
   Rot[8]=(RotPlane==YZ||RotPlane==XZ?CosA:1);
   switch(RotPlane){
      case(XY):{
         Rot[1]=-SinA;
         Rot[3]=SinA;
         break;
      }
      case(XZ):{
         Rot[2]=SinA;
         Rot[6]=-SinA;
         break;
      }
      case(YZ):{
         Rot[5]=-SinA;
         Rot[7]=SinA;
         break;
      }
   }
}

void GeomManipulator::Rotate(const double Theta,const Plane RotPlan){
   double Rot[9];
   memset(&Rot[0],0.0,sizeof(double)*9);
   RotMat(Rot,Theta,RotPlan);
   Rotate(Rot);
}

void GeomManipulator::Set(){
   SetCarts(Carts_);
}
}}//End namespaces
