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
#include <set>
#include <math.h>
#include "UnitCellGuts.h"
#include "Units.h"
#include "../Utils/GeomManipulator.h"
#include "../Utils/Geometry.h"
namespace psi {
namespace LibMolecule {
typedef std::vector<int>::iterator IntIt;
inline void GenericConverter(double ConvFac, double param[], const double p1,
      const double p2, const double p3) {
   param[0]=p1*ConvFac;
   param[1]=p2*ConvFac;
   param[2]=p3*ConvFac;
}

void UnitCellGuts::SetAngles(const double alpha, const double beta,
      const double gamma, const bool IsDegree) {
   BaseUnitConverter Conv;
   GenericConverter(IsDegree ? Conv(DEGREE, RADIAN) : 1.0, angles_, alpha, beta,
         gamma);

}

void UnitCellGuts::SetSides(const double a, const double b, const double c,
      const bool IsBohr) {
   BaseUnitConverter Conv;
   GenericConverter(IsBohr ? 1.0 : Conv(ANGSTROM, BOHR), sides_, a, b, c);
}

/*
 void Shift(const double* TVec, GeomManipulator& Manip, std::vector<int>& NoBond,
 Molecule* Mol) {
 Geometry NewGeom(Mol);
 std::vector<int> NewNoBond;
 for (IntIt AtomI=NoBond.begin(); AtomI!=NoBond.end(); ++AtomI) {
 if (NewGeom.NBonds((*AtomI))==0) {
 NewNoBond.push_back((*AtomI));
 if (TVec!=NULL) Manip.Translate(&TVec[0], (*AtomI));
 }
 }
 Manip.Set();
 NoBond=NewNoBond;
 }
 *
 *
 * bool UnitCellGuts::FixUnitCell(std::vector<int>& NoBond, int Direction) {
 GeomManipulator Manip(this);
 std::vector<double> Translate(3, 0.0),TVec(3, 0.0);
 if (Direction!=-1) {
 Translate[Direction]=1.0;
 DGEMM_Int(&TVec[0], Frac2Cart_, &Translate[0], 3, 3, 1);
 }
 Shift((Direction==-1 ? NULL : &TVec[0]), Manip, NoBond, this);
 if (NoBond.size()>0&&Direction!=-1) {
 //Minus 1 moves it back, Minus 2 tries the opposite direction
 for (int i=0; i<3; i++)
 TVec[i]*=-2.0;
 Shift(&TVec[0], Manip, NoBond, this);
 if (NoBond.size()>0) {
 //Still didn't work move it back to the cell try next direction
 for (int i=0; i<3; i++)
 TVec[i]*=-0.5;
 Shift(&TVec[0], Manip, NoBond, this);
 }
 }
 return NoBond.size()==0;
 }*/



UnitCellGuts::UnitCellGuts(const Molecule& Mol, const double* Sides,
      const double* Angles, const bool IsFrac, const bool IsBohr,
      const bool IsDegree) :
      Molecule(Mol) {
   SetSides(Sides[0], Sides[1], Sides[2], IsBohr);
   SetAngles(Angles[0], Angles[1], Angles[2], IsDegree);
   SetTrans();
   if (IsFrac) {
      GeomManipulator Manip(this);
      //Unscale the fractional coordinates since the conversion wasn't
      //appropriate
      Manip.Scale(1/AngToBohr());
      Manip.Rotate(Frac2Cart_);
      Manip.Set();
   }
}

void UnitCellGuts::Copy(const UnitCellGuts& other) {
   for (int i=0; i<3; i++) {
      this->angles_[i]=other.angles_[i];
      this->sides_[i]=other.sides_[i];
      for (int j=0; j<3; j++) {
         this->Frac2Cart_[i*3+j]=other.Frac2Cart_[i*3+j];
         this->Cart2Frac_[i*3+j]=other.Cart2Frac_[i*3+j];
      }
   }
}

//-----The def below is just the two transformation matrices hard-coded----
//It's stolen from Wikipedia
void UnitCellGuts::SetTrans() {
   double costheta[3],sintheta[3],cos2=0.0,cosprod=1.0;
   for (int i=0; i<3; i++) {
      costheta[i]=cos(angles_[i]);
      cos2+=costheta[i]*costheta[i];
      cosprod*=costheta[i];
      sintheta[i]=sin(angles_[i]);
   }
   double volume=sqrt(1-cos2+2*cosprod);
//Easy elements
   Frac2Cart_[0]=sides_[0];
   Cart2Frac_[0]=1.0/sides_[0];
   Frac2Cart_[3]=Cart2Frac_[3]=0.0;
   Frac2Cart_[6]=Cart2Frac_[6]=0.0;
   Frac2Cart_[7]=Cart2Frac_[7]=0.0;

//More difficult elements
   double cosbg=costheta[1]*costheta[2],invsin2=1/sintheta[2];
   Frac2Cart_[1]=sides_[1]*costheta[2]; //b*cos(gamma)
//-cos(gamma)/a*sin(gamma)
   Cart2Frac_[1]=-costheta[2]*Cart2Frac_[0]*invsin2;
   Frac2Cart_[2]=sides_[2]*costheta[1];   //c*cos(beta)
//cos(alpha)*cos(gamma)-cos(beta)/av*sin(gamma)
   Cart2Frac_[2]=Cart2Frac_[0]*invsin2*(costheta[0]*costheta[2]-costheta[1])/volume;
   Frac2Cart_[4]=sides_[1]*sintheta[2];   // b*sin(gamma)
   Cart2Frac_[4]=1/Frac2Cart_[4];   //1/b*sin(gamma)
//c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
   Frac2Cart_[5]=sides_[2]*invsin2*(costheta[0]-cosbg);
//cos(beta)cos(gamma)-cos(alpha)/bv*sin(gamma)
   Cart2Frac_[5]=invsin2*(cosbg-costheta[0])/(sides_[1]*volume);
   Frac2Cart_[8]=sides_[2]*volume*invsin2;   // cv/sin(gamma)
   Cart2Frac_[8]=1/Frac2Cart_[8];   //sin(gamm)/cv
}

}}//End namespaces

