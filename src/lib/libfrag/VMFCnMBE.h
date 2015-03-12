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
#ifndef SRC_LIB_LIBFRAG_VMFCNMBE_H_
#define SRC_LIB_LIBFRAG_VMFCNMBE_H_

#include "CanonicalMBE.h"
#include "libmolecule/Utils/BSSEer.h"
#include "libPsiUtil/Table.h"
namespace psi {
namespace LibFrag {

/** \brief Handles the MBE, when VMFCn BSSE corrections are applied
 *
 *  When a VMFCn BSSE correction is applied the canonical MBE changes
 *  form in that we can't use monomer energies to compute two-body
 *  corrections anymore (we need monomers in the dimer basis set).
 *  Similar ideas hold for higher orders.  Note that we have all
 *  the information necessary to compute the canonical MBE, and so
 *  this equation does that as well.
 *
 */
class VMFCnMBE:public CanonicalMBE {
   private:
      template <typename T>
      void VMFCnEgy(const int N,
                    const LibMolecule::FragmentedSystem& Systems,
                    const MBEProp<T>&MonoProperties,
                    MBEProp<T>& CanonicalEgys,
                    MBEProp<T>& VMFCEgys)const;
      /*template <typename T>
      void ApproxVMFCnEgy(const int N,
                    const LibMolecule::FragmentedSystem& Systems,
                    const MBEProp<T>&MonoProperties,
                    MBEProp<T>& ApproxEgys)const;*/
   public:
      VMFCnMBE(int N=1) :
            CanonicalMBE(N) {
      }
      ///Computes and returns the property
      template <typename T>
      MBEProp<T> PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
            const MBEProp<T>& MonoProperties) ;
};

/***************Implementations are below********************/

template <typename T>
std::string PrintOut(
              const std::string& Name,
              const MBEProp<T>& CanEgys,
              const MBEProp<T>& VMFCEgys);

template <typename T>
MBEProp<T> VMFCnMBE::PropertyImpl(const LibMolecule::FragmentedSystem& Systems,
      const MBEProp<T>&MonoProperties)  {
   //Start by computing the normal MBE energies
   MBEProp<T> CanonicalEgys=CanonicalMBE::PropertyImpl(Systems, MonoProperties);
   MBEProp<T> VMFCEgys(N, "VMFC(n) Energy");
   //MBEProp<T> ApproxEgys(N,"Approx VMFC(n)");
   VMFCnEgy(N,Systems,MonoProperties,CanonicalEgys,VMFCEgys);
   //ApproxVMFCnEgy(N,Systems,MonoProperties,ApproxEgys);
   LastResult_=PrintOut(MonoProperties.Name(),CanonicalEgys,VMFCEgys);
   return VMFCEgys;
}

/** \brief Computes the n-body interaction of a fragment
 *
 *   We assume we were given:
 *   \verbatim
 *          (M+1)...n
 *   IJK...M
 *   \endverbatim
 *
 *   The interaction of this is just it's energy less the interaction of
 *   each (M-1)-mer in the n-mer basis, i.e.
 *   \verbatim
 *              M...n
 *   IJK...(M-1)
 *   \endverbatim
 *   down to:
 *   \verbatim
 *    (M+1)...n
 *   M
 *   \endverbatim
 *
 *
 *
 */
template<typename T>
void Interaction(const LibMolecule::SerialNumber& SN,
      const MBEProp<T>& Props,MBEProp<T>& Int){
   typedef LibMolecule::SerialNumber SN_t;
   SN_t FakeSN,Real;FakeSN.insert(0);
   Int.Change(0,FakeSN,Props(SN.size()-1,SN));
   SN_t::const_iterator SNJ=SN.begin(),SNJEnd=SN.end();
   for(;SNJ!=SNJEnd;++SNJ)
      if((*SNJ)>0)Real.insert(*SNJ);
   if(Real.size()>1){
      LibMolecule::VMFCnItr SNI(SN,false);
      for(;!SNI.Done();++SNI){
         /*SN_t Real2;
         SNJ=SNI->begin(),SNJEnd=SNI->end();
         for(;SNJ!=SNJEnd;++SNJ)
            if((*SNJ)>0)Real2.insert((*SNJ));
         double sign=((Real.size()+Real2.size())%2==1?-1.0:1.0);*/
         MBEProp<T> Result(1);
         Interaction((*SNI),Props,Result);
         Int.Change(0,FakeSN,Result(0,FakeSN)*-1.0);
      }
   }
}

template <typename T>
void VMFCnMBE::VMFCnEgy(const int N,const LibMolecule::FragmentedSystem& Systems,
              const MBEProp<T>&MonoProperties,
              MBEProp<T>& CanonicalEgys,
              MBEProp<T>& VMFCEgys)const{
   typedef LibMolecule::FragmentedSystem Frags_t;
   typedef LibMolecule::SerialNumber SN_t;
   SN_t FakeSN;
   FakeSN.insert(0);
   VMFCEgys.Change(0, FakeSN, CanonicalEgys(0, FakeSN));
   for (int n=1; n<N; n++) {
      VMFCEgys.Change(n, FakeSN, VMFCEgys(n-1, FakeSN));
      Frags_t::iterator FragI=Systems.begin(n),FragIEnd=Systems.end(n);
      for (; FragI!=FragIEnd; ++FragI) {
         if (!IsReal(FragI)) break;
         SN_t SN=(*FragI)->GetSN();
         MBEProp<T> Result(1);
         Interaction(SN,MonoProperties,Result);
         VMFCEgys.Change(n,FakeSN,Result(0,FakeSN));
      }
   }
}

/*
template <typename T>
void VMFCnMBE::ApproxVMFCnEgy(const int N,
              const LibMolecule::FragmentedSystem& Systems,
              const MBEProp<T>&MonoProperties,
              MBEProp<T>& ApproxEgys)const{
   typedef LibMolecule::SerialNumber SN_t;
   SN_t FakeSN;FakeSN.insert(0);
   LibMolecule::FragmentedSystem::iterator NMerI=Systems.begin(N-1),
                                        NMerIEnd=Systems.end(N-1);
   MBEProp<T> Interactions(N);
   PsiMap<SN_t,unsigned int> Counter;
   for(;NMerI!=NMerIEnd;++NMerI){
      SN_t  SN=(*NMerI)->GetSN();
      if(*SN.begin()<0)break;
      MBEProp<T> Result(1);
      Interaction(SN,MonoProperties,Result);
      Interactions.Change(N-1,FakeSN,Result(0,FakeSN));
      LibMolecule::VMFCnItr SNI(SN,false);
      for(;!SNI.Done();++SNI){
         SN_t Real;
         SN_t::const_iterator FragI=SNI->begin(),FragIEnd=SNI->end();
         for(;FragI!=FragIEnd;++FragI)
            if((*FragI)>0)Real.insert((*FragI));
         if(Counter.count(Real)==0)Counter[Real]=1;
         else ++Counter[Real];
         MBEProp<T> Result2(1);
         Interaction((*SNI),MonoProperties,Result2);
         Interactions.Change(Real.size()-1,Real,Result2(0,FakeSN));
      }
   }
   for(int i=0;i<N;i++){
      typename MBEProp<T>::const_iterator IntJ=Interactions.begin(i),
            IntJEnd=Interactions.end(i);
      if(i>0)ApproxEgys.Change(i,FakeSN,ApproxEgys(i-1,FakeSN));
      for(;IntJ!=IntJEnd;++IntJ){
            ApproxEgys.Change(i,FakeSN,IntJ->second*
                  (i<N-1?1/Counter[IntJ->first]:1.0));
      }
   }
}
*/


static std::string PrintDouble(double x){
   std::stringstream Temp1;
   Temp1<<std::setw(14)<<std::setprecision(14)<<std::fixed<<x;
   return Temp1.str();
}

template <typename T>
std::string PrintOut(const std::string& PropName,
              const MBEProp<T>& CanEgys,
              const MBEProp<T>& VMFCEgys){
   std::vector<std::string> Values;
   LibMolecule::SerialNumber  FakeSN;FakeSN.insert(0);
   const int NValues=4;
   for(int i=0;i<CanEgys.N();i++){
      std::stringstream Temp,Temp1,Temp2,Temp3;
      Temp<<i+1<<"-body";
      Values.push_back(Temp.str());
      Values.push_back(PrintDouble(CanEgys(i,FakeSN)));
      Values.push_back(PrintDouble(VMFCEgys(i,FakeSN)));
      Values.push_back(PrintDouble(CanEgys(i,FakeSN)-VMFCEgys(i,FakeSN)));
   }
   typedef TableColumn<std::string> StrCol;
   std::vector<StrCol> Columns;
   Columns.push_back(StrCol(PropName,&Values[0],NValues,0,'\0','|'));
   Columns.push_back(StrCol("Canonical E",&Values[1],NValues));
   Columns.push_back(StrCol("VMFC(n) E",&Values[2],NValues));
   Columns.push_back(StrCol("BSSE",&Values[3],NValues));
   Table<StrCol,StrCol,StrCol,StrCol> MyTable(
         CanEgys.N(),Columns[0],Columns[1],Columns[2],Columns[3]);
   MyTable.SetBorder(TOP,'*');
   MyTable.SetBorder(BOTTOM,'*');
   return MyTable.GetTable();
}


}}      //End namespaces

#endif /* SRC_LIB_LIBFRAG_VMFCNMBE_H_ */
