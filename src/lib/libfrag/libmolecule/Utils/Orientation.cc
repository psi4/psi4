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

#include "Orientation.h"

namespace psi{
namespace LibMolecule{

void Orientation::Reorient(Molecule& Mol){
   GeomManipulator Manip(&Mol);
   CalcOrigin(Mol);
   Manip.Translate(Trans2Origin_,-1.0);
   SetZAxis(Mol,Manip);
   Manip.Rotate(Rotation_);
   double OrientRot[]={1.0,0.0,0.0,
                       0.0,1.0,0.0,
                       0.0,0.0,1.0};
   SetDirection(Mol,Manip,OrientRot);
   Manip.Rotate(OrientRot);
   double Theta=GetTheta(Mol,Manip);
   if(Theta!=0.0)Manip.Rotate(Theta,XY);
   Manip.Set();
}

///Function that computes the center of value-weighted center
void EigenOrientation::CalcOrigin(const Molecule& other) {
   double TotalValue=0.0;
   MolItr AtomI=other.Begin(),AtomEnd=other.End();
   for (; AtomI!=AtomEnd; ++AtomI) {
      const double ValueI=Value((*AtomI).get());
      TotalValue+=ValueI;
      for (int i=0; i<3; i++)Trans2Origin_[i]+=(*(*AtomI))[i]*ValueI;
   }
   for (int i=0; i<3; i++)Trans2Origin_[i]/=TotalValue;
}


void EigenOrientation::SetZAxis(const Molecule& other,
                 const GeomManipulator& Manip) {
   MolItr AtomI=other.Begin(),AtomEnd=other.End();
   for (int i=0; AtomI!=AtomEnd; ++AtomI, ++i) {
      const double Z=Value((*AtomI).get()),
                  x2=Manip(i, 0)*Manip(i, 0),
                  y2=Manip(i, 1)*Manip(i, 1),
                  z2=Manip(i, 2)*Manip(i, 2),
                  Zx=Z*Manip(i, 0);
      Rotation_[0]+=Z*(y2+z2);
      Rotation_[1]-=Zx*Manip(i, 1);
      Rotation_[2]-=Zx*Manip(i, 2);
      Rotation_[4]+=Z*(x2+z2);
      Rotation_[5]-=Z*Manip(i, 1)*Manip(i, 2);
      Rotation_[8]+=Z*(x2+y2);
   }
   Rotation_[3]=Rotation_[1];
   Rotation_[6]=Rotation_[2];
   Rotation_[7]=Rotation_[5];
   //MoC'=U[MoC]U^dagger and MoC=U on return
   Diagonalize(Rotation_, 3, Omega_);
   //Eigenvalues are in ascending order and necessarily positive
   const double d1=Omega_[1]-Omega_[0],
                d2=Omega_[2]-Omega_[1];
   const bool D1Degen=(d1<Tolerance_),
              D2Degen=(d2<Tolerance_),
              linear=(Omega_[0]<Tolerance_||
                      Omega_[1]<Tolerance_||
                      Omega_[2]<Tolerance_);
   //For symmetric tops, linear molecules the Z-axis is the unique one
   //and showed up as the X-Axis under the following two conditions:
   if ((!D1Degen&&D2Degen)||linear) {
      //Need to switch eigenvector 1 to eigenvector 3
      for (int i=0; i<3; i++) {
         double temp=Rotation_[0*3+i];
         Rotation_[0*3+i]=Rotation_[2*3+i];
         Rotation_[2*3+i]=temp;
      }
   }

}

void EigenOrientation::SetDirection(const Molecule& Mol,
      const GeomManipulator& Manip,double *Rot){
   double FourthMoments[]={0.0,0.0,0.0,
                           0.0,0.0,0.0};
   MolItr AtomI=Mol.Begin(),AtomEnd=Mol.End();
   for(int i=0;AtomI!=AtomEnd;++AtomI,++i){
      const double Z=Value((*AtomI).get());
      for(int j=0;j<3;j++){
         const double q2=Manip(i,j)*Manip(i,j);
         FourthMoments[(Manip(i,j)>=0.0?j:j+3)]+=Z*q2*q2;
      }
   }
   for(int i=0;i<3;i++)
      if(FourthMoments[i+3]<FourthMoments[i]){
         Rot[i*3+i]=-1.0;
         for(int j=0;j<3;j++)Rotation_[j*3+i]*=-1.0;
      }
}

double Distance2Z(const GeomManipulator& Manip, const int I){
   double d=0.0;
   for(int i=0;i<2;i++)d+=Manip(I,i)*Manip(I,i);
   return sqrt(d);
}

void UpdateParams(int i,std::vector<int>& FoundAtoms,
                  double& SmallestAbsZ,bool& pos,
                  double& distance2Z,int& Zs,
                  const MolItr& AtomI,const GeomManipulator& Manip){
   FoundAtoms.clear();
   FoundAtoms.push_back(i);
   SmallestAbsZ=fabs(Manip(i, 2));
   pos=(Manip(i,2)>0);
   distance2Z=Distance2Z(Manip,i);
   Zs=(*AtomI)->Z();
}

std::vector<int> SphericalTop(const GeomManipulator& Manip,
      const Molecule& other,
      const double Tolerance){
   double ShortestD=10000;
   std::vector<int> AtomI;
   for (int i=0; i<other.NAtoms(); i++) {
      double d=0.0;
      for (int j=0; j<3; j++)d+=Manip(i, j)*Manip(i, j);
      d=sqrt(d);
      double dz=ShortestD-d;
      ///Don't let the origin dictate this
      if (dz>Tolerance&&d>Tolerance){
         ShortestD=d;
         AtomI.clear();
         AtomI.push_back(i);
      }
      else if(fabs(dz)<Tolerance&&d>Tolerance){
         AtomI.push_back(i);
      }
   }
   return AtomI;
}


std::vector<int> SymmetricTop(const GeomManipulator& Manip,
                             const Molecule& other,
                              const double Tolerance) {
   std::vector<int> FoundAtoms;
   int Zs=-1;
   bool pos=true;
   double SmallestAbsZ=100000,distance2Z=10000;
   MolItr AtomI=other.Begin(),AtomEnd=other.End();
   for (int i=0;AtomI!=AtomEnd; i++,++AtomI) {
      //Skip atoms on the Z-Axis
      if (fabs(Manip(i, 0))<Tolerance&&fabs(Manip(i, 1))<Tolerance) continue;
      double dZ=SmallestAbsZ-fabs(Manip(i, 2));
      if (dZ>Tolerance)
         UpdateParams(i,FoundAtoms,SmallestAbsZ,pos,distance2Z,Zs,AtomI,Manip);
      ///Fabs prevents the cases where Manip(i,2)>SmallestAbsZ
      else if (fabs(dZ)<=Tolerance){
         //If we've been finding negative z's, but now found a positive
         //one, we use it instead
         if(!pos&&Manip(i,2)>0)
            UpdateParams(i,FoundAtoms,SmallestAbsZ,pos,distance2Z,Zs,AtomI,Manip);
         //Similarly, if we've been finding positive z's and now found
         //a negative ignore it
         else if(pos&&Manip(i,2)<0)continue;
         //Found an atom with the same sign, and Z component
         else{
            double tempdist=Distance2Z(Manip,i);
            double dz=tempdist-distance2Z;
            //Is the new atom closer?
            if(fabs(dz)>Tolerance&&dz<0.0){//It is
               UpdateParams(i,FoundAtoms,SmallestAbsZ,pos,distance2Z,Zs,AtomI,Manip);
            }
            else if(fabs(dz)>Tolerance&&dz>0.0)continue;//Further
            else{
               //Does the new atom have a lower Z
               if((*AtomI)->Z()<Zs)
                  UpdateParams(i,FoundAtoms,SmallestAbsZ,
                        pos,distance2Z,Zs,AtomI,Manip);
               else FoundAtoms.push_back(i);//Same circular set
            }

         }
      }
   }
   return FoundAtoms;
}

double EigenOrientation::GetTheta(const Molecule& Mol,
                                  const GeomManipulator& Manip)const{
   const double d1=Omega_[1]-Omega_[0],
                d2=Omega_[2]-Omega_[1];
   const bool D1Degen=(d1<Tolerance_),
              D2Degen=(d2<Tolerance_),
              linear=(Omega_[0]<Tolerance_||
                      Omega_[1]<Tolerance_||
                      Omega_[2]<Tolerance_);
   //At this point if omega[0]!=omega[1]!=omega[2], we are done
   //Also done for linear molecules, and atoms.
   if ((D1Degen||D2Degen)&&!linear) {
      bool IsSpherical=D1Degen&&D2Degen;
      std::vector<int> FoundAtoms=(IsSpherical?
            SphericalTop(Manip,Mol,Tolerance_):
            SymmetricTop(Manip,Mol,Tolerance_));
      double YAxis[]={0.0,1.0,0.0};
      double dx=Manip(FoundAtoms[0],0),dy=Manip(FoundAtoms[0],1);
      double CosTheta=dy/sqrt(dx*dx+dy*dy);
      return acos(CosTheta);
   }
   return 0.0;
}

double MassOrientation::Value(const Atom* I)const{
   return I->Mass();
}

double ChargeOrientation::Value(const Atom* I)const{
   return I->Z();
}
}}//End namespaces


