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

#ifndef GMBE_H_
#define GMBE_H_
#include "FragmentedSys.h"
#include "MBEProp.h"
namespace psi{
namespace LibFrag{

template<typename T>
class Expansion{
   private:
      T Impl;
   public:
      Expansion(const int N):Impl(N){}
      template<typename T2>
      MBEProp<T2> Property(const LibMolecule::FragmentedSystem& Systems,
            const MBEProp<T2>& MonoProperties){
         return Impl.PropertyImpl(Systems,MonoProperties);
      }
      std::string PrintOut()const{return Impl.PrintOutImpl();}

};

class ExpanImplBase{
   protected:
      int N;
      double Phase(const int i)const{return (i%2==0?1:-1);}
      ///Once calculated this will be the result in a pretty format
      std::string LastResult_;
   public:
      ExpanImplBase(const int NewN):N(NewN){}
      ///Returns True, useful if for some reason you want to know this
      virtual bool IsGMBE()const{return true;}
      virtual ~ExpanImplBase(){}
      ///Returns the result in a pretty format
      std::string PrintOutImpl()const{return LastResult_;}
};

class GMBE:public ExpanImplBase{
	public:
		///Makes an N-body GMBE, N defaults to 1
		GMBE(int newN=1):ExpanImplBase(newN){}

		///Returns the total energy of this GMBE expansion
		template<typename T>
		MBEProp<T> PropertyImpl(
		      const LibMolecule::FragmentedSystem& Systems,
		      const MBEProp<T>& MonoProperties);
};

template<typename T>
MBEProp<T> GMBE::PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
      const MBEProp<T>& MonoProperties){
   MBEProp<T> Prop(1);
   /*If N==1 special case and energies are only in Energies[0],but
   //we get sizes from Systems[1] and Systems[2]
   int index1=(N==1?0:1);
   for (int i=0; i<Systems[1].size(); i++)
      TEnergy+=Systemsts[i]*Energies[index1][i];

   for (int j=0; j<Systems[2].size(); j++) {
      int index2=(N==1?Systems[1].size()+j:j);
      TEnergy-=NMults[j]*Energies[index1][index2];
   }
   psi::outfile->Printf( "The total %d-body GMBE %s is: %16.12f (a.u.)", N,
         RealName.c_str(),TEnergy);*/
   return Prop;
}


}}//End namespaces



#endif /* GMBE_H_ */
