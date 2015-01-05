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

#include "UnitCell.h"
#include "Units.h"
namespace psi{
namespace LibMolecule{

inline void GenericConverter(double ConvFac,double param[],const double p1,
                      const double p2, const double p3){
   param[0]=p1*ConvFac;
   param[1]=p2*ConvFac;
   param[2]=p3*ConvFac;
}

void UnitCell::SetAngles(const double alpha, const double beta,
                         const double gamma, const bool IsDegree){
   BaseUnitConverter Conv;
   GenericConverter(IsDegree?Conv(DEGREE,RADIAN):1.0,angles_,
                    alpha,beta,gamma);

}

void UnitCell::SetSides(const double a, const double b, const double c,
                        const bool IsBohr){
   BaseUnitConverter Conv;
   GenericConverter(IsBohr?1.0:Conv(ANGSTROM,BOHR),sides_,a,b,c);
}


//-----The def below is just the two transformation matrices hard-coded----
void UnitCell::SetTrans(){
   double cos[3],sin[3],cos2=0.0,cosprod=1.0;
   for(int i=0;i<3;i++){
      cos[i]=std::cos(angles_[i]);
      cos2+=cos[i]*cos[i];
      cosprod*=cos[i];
      sin[i]=std::sin(angles_[i]);
   }
   double volume=sqrt(1-cos2+2*cosprod);
   //Easy elements
   Frac2Cart_[0]=sides_[0];Cart2Frac_[0]=1.0/sides_[0];
   Frac2Cart_[3]=Cart2Frac_[3]=0.0;
   Frac2Cart_[6]=Cart2Frac_[0]=0.0;
   Frac2Cart_[7]=Cart2Frac_[0]=0.0;

   //More difficult elements
   double cosbg=cos[1]*cos[2],
           ctimesvol=sides_[2]*volume,
           invsin2=1/sin[2];
   Frac2Cart_[1]=sides_[1]*cos[2];//b*cos(gamma)
   //-cos(gamma)/a*sin(gamma)
   Cart2Frac_[1]=-cos[2]*Cart2Frac_[0]*invsin2;
   Frac2Cart_[2]=sides_[2]*cos[1];//c*cos(beta)
   //cos(alpha)*cos(gamma)-cos(beta)/av*sin(gamma)
   Cart2Frac_[2]=Cart2Frac_[0]*invsin2*(cos[0]*cos[2]-cos[1])/volume;
   Frac2Cart_[4]=sides_[1]*sin[2];// b*sin(gamma)
   Cart2Frac_[4]=invsin2/sides_[1];//1/b*sin(gamma)
   //c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
   Frac2Cart_[5]=sides_[2]*invsin2*(cos[0]-cosbg);
   //cos(beta)cos(gamma)-cos(alpha)/bv*sin(gamma)
   Cart2Frac_[6]=invsin2*(cosbg-cos[0])/(sides_[1]*volume);
   Frac2Cart_[8]=ctimesvol*invsin2;// cv/sin(gamma)
   Cart2Frac_[8]=sin[2]/ctimesvol;//sin(gamm)/cv
}


}}

