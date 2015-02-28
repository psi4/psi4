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

#ifndef MBE_H_
#define MBE_H_
#include <set>
#include <utility>
#include "GMBE.h"
#include <boost/math/special_functions/binomial.hpp>
namespace psi {
namespace LibFrag {

/** \brief Specialization of the GMBE to non-intersecting fragments
 *
 *  For the case of the MBE an \f$n\f$-body property \f$P\f$ of an
 *  \f$N\f$ fragment system can be written:
 *  \f[
 *  P\approx\sum_{i=1}^n(-1)^{n-i}{_{N-i-1}}C_{n-i}P^{(i)},
 *  \f]
 *  where \f${_{\alpha}}C_{\beta}\f$ is \f$\alpha\f$ choose \f$\beta\f$
 *  and \f$P^{(i)}\f$ is the sum of the property over all $i$-mers.
 *
 */
class MBE:public ExpanImplBase {
   private:
      ///Returns N choose M
      double Binomial(const int N, const int M)const {
         return boost::math::binomial_coefficient<double>(N, M);
      }
      ///Returns the total N-body property, assuming nfrags and m<N
      template <typename T>
      MBEProp<T> NBodyProp(const int N, const int nfrags,
            const MBEProp<T>& DeltaEs)const;
      ///Returns the coefficient for the m-th term of an N-body expansion of nfrags.  m=1 for Ens[0], m=N for Ens[N-1]
      double Coef(const int nfrags, const int N, const int m) const{
         return Phase(N-m)*Binomial(nfrags-m-1, N-m);
      }
   public:
      ///Sets MBE truncation order to N, default is 1
      MBE(int N=1) :ExpanImplBase(N) {
      }
      ///Computes and returns the energy
      template <typename T>
      MBEProp<T> PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
            const MBEProp<T>& MonoProperties)const;
      bool IsGMBE() const {return false;}
};

///Does the MBE by actually calculating the \f$\Delta E_{ij\cdots n}\f$ terms
class CanonicalMBE:public ExpanImplBase{
   private:
      std::vector<LibMolecule::SerialNumber::const_iterator> FillIndices(
            const int m,
            const LibMolecule::SerialNumber& CurrSerial)const;
      bool UpdateIndices(
            const LibMolecule::SerialNumber::const_iterator& LastIndex,
            const int m,
            std::vector<LibMolecule::SerialNumber::iterator>& Indices)const;
   public:
      CanonicalMBE(int N=1):ExpanImplBase(N) {}
      ///Computes and returns the property
      template<typename T>
      MBEProp<T> PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
            const MBEProp<T>& MonoProperties)const;
};


template <typename T>
MBEProp<T> MBE::PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
      const MBEProp<T>& MonoProperties)const {
   typedef typename MBEProp<T>::const_iterator Itr_t;
   //This is the value of the property at all levels up to including N
   MBEProp<T> Value(N);
   //These will be the corrections to the property
   MBEProp<T> Corrs(N);
   LibMolecule::SerialNumber Dummy;
   Dummy.insert(0);
   int TrueNumMono=0;
   LibMolecule::FragmentedSystem::iterator MonoI=Systems.begin(0),
                                           MonoEnd=Systems.end(0);
   for(;MonoI!=MonoEnd;++MonoI)
      TrueNumMono+=Systems.Coef(0,(*MonoI)->GetSN());
   //Total property of each set of n-mers
   MBEProp<T> En(N);
   for (int i=0; i<N; i++) {
      Itr_t PropI=MonoProperties.begin(i),PropEnd=MonoProperties.end(i);
      MonoI=Systems.begin(i);
      for (; PropI!=PropEnd; ++PropI)
         En.Change(i, Dummy, (*PropI)->second*Systems.Coef(i,(*MonoI)->GetSN()));
      Value.Change(i, Dummy, NBodyProp(i+1, TrueNumMono, En)(0, 0));
      if (i==0) continue;
      Corrs.Change(i, Dummy, Value(i,0));
      Corrs.Change(i, Dummy, Value(i-1,0)*-1.0);
   }
   return Value;
}

template <typename T>
MBEProp<T> MBE::NBodyProp(const int N, const int nfrags,
      const MBEProp<T>& DeltaEs)const {
   LibMolecule::SerialNumber Dummy;
   Dummy.insert(0);
   MBEProp<T> Result(1,"Total Energy");
   if (N!=nfrags){
      for (int i=0; i<N; i++){
         Result.Change(0, Dummy, DeltaEs(i, Dummy)*MBE::Coef(nfrags, N, i+1));
      }
   }
   else Result.Change(0,Dummy,DeltaEs(N-1,Dummy));
   return Result;
}


template<typename T>
MBEProp<T> CanonicalMBE::PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
      const MBEProp<T>& MonoProperties)const{
   typedef LibMolecule::SerialNumber Set_t;
   typedef typename MBEProp<T>::const_iterator Itr_t;
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
         const Set_t& CurrSerial=(*NMerI)->GetSN();
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
   return Value;
}


}}//End namespaces

#endif /* MBE_H_ */
