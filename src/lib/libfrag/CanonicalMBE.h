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
#ifndef SRC_LIB_LIBFRAG_CANONICALMBE_H_
#define SRC_LIB_LIBFRAG_CANONICALMBE_H_
#include "MBE.h"

namespace psi{
namespace LibFrag{

///Does the MBE by actually calculating the \f$\Delta E_{ij\cdots n}\f$ terms
class CanonicalMBE:public ExpanImplBase{
   protected:
      virtual std::vector<LibMolecule::SerialNumber::const_iterator> FillIndices(
            const int m,
            const LibMolecule::SerialNumber& CurrSerial)const;
      virtual bool UpdateIndices(
            const LibMolecule::SerialNumber::const_iterator& LastIndex,
            const int m,
            std::vector<LibMolecule::SerialNumber::iterator>& Indices)const;
      bool IsReal(const LibMolecule::FragmentedSystem::iterator& Itr)const;
   public:
      CanonicalMBE(int N=1):ExpanImplBase(N) {}
      ///Computes and returns the property
      template<typename T>
      MBEProp<T> PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
            const MBEProp<T>& MonoProperties);
};

template<typename T>
MBEProp<T> CanonicalMBE::PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
      const MBEProp<T>& MonoProperties){
   typedef LibMolecule::SerialNumber Set_t;
   Set_t Dummy;
   Dummy.insert(0);
   //This is the value of the property at all levels up to including N
   MBEProp<T> Value(N,MonoProperties.Name());
   //These will be the corrections to the property
   MBEProp<T> Corrs(N), Properties(N);
   for(int n=0;n<N;n++){
      if(n>0)Value.Change(n,Dummy,Value(n-1,Dummy));
      LibMolecule::FragmentedSystem::iterator NMerI=Systems.begin(n),
                                             NMerIEnd=Systems.end(n);
      for(;NMerI!=NMerIEnd;++NMerI){
         if(!IsReal(NMerI))break;
         const LibMolecule::SerialNumber& CurrSerial=(*NMerI)->GetSN();
         //Last index in the set
         Set_t::const_iterator LastIndex=CurrSerial.end();
         MBEProp<T> Prop(1);
         Prop.Change(0,Dummy,MonoProperties[n][CurrSerial]);
         for(int m=n-1;m>=0;m--){
            std::vector<Set_t::const_iterator> Indices=
                                                   FillIndices(m,CurrSerial);
            bool done=false;
            do{
               //Read the indices into a set
               Set_t TempSet,TempTempSet;
               for(int l=0;l<=m;l++)TempTempSet.insert((*Indices[l]));
               TempSet=Systems.SNLookUp(TempTempSet);
               //Now we found our TempSet
               Prop.Change(0,Dummy,Properties[m][TempSet]*-1.0);
               bool goodindex=UpdateIndices(LastIndex,m,Indices);
               if(!goodindex)done=true;
            }while(!done);
         }
         const Set_t& SN=Systems.SNLookUp(CurrSerial);
         Set_t::const_iterator SetI=CurrSerial.begin(),
                               SetEnd=CurrSerial.end();
         Properties[n][CurrSerial]=Prop(0,Dummy);
         if(Systems.Coef(n,SN)!=0.0){
            Corrs.Change(n,Dummy,Prop(0,Dummy)*Systems.Coef(n, SN));
            Value.Change(n,Dummy,Prop(0,Dummy)*Systems.Coef(n, SN));
         }
      }
   }
   LastResult_=Value.PrintOut();
   return Value;
}

}}//End namespaces


#endif /* SRC_LIB_LIBFRAG_CANONICALMBE_H_ */
