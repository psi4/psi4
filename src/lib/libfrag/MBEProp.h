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
#ifndef SRC_LIB_LIBFRAG_MBEPROP_H_
#define SRC_LIB_LIBFRAG_MBEPROP_H_
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "libmolecule/LibFragFragment.h"
#include "libPsiUtil/PsiMap.h"
#include "exception.h"
namespace psi{
namespace LibFrag{

/** \brief A class to hold generic Many-Body properties
 *
 *   The template parameter T is the type of the property and that
 *   property defines:
 *   \code
 *   const T& operator+=(const T& other);
 *   const T& operator*(const double scalar)const;
 *   \endcode
 *
 *   then the MBE and GMBE will be able to reassemble your property,
 *   automatically for you. The PrintOut fxn also currently assumes
 *   that you can pass your property to a stringstream (possible by
 *   overloading the ostream::operator<< operator for your type);
 *   however, I'm not wed to this and this fxn can be commented out
 *   for more complex template parameters
 */

template<typename T>
class MBEProp:protected std::vector<PsiMap<LibMolecule::SerialNumber,T> >{
   private:
      typedef PsiMap<LibMolecule::SerialNumber,T> BaseArg_t;
      ///Convenient definition of the base type
      typedef std::vector<BaseArg_t> Base_t;
      ///What this property is called (energy, gradient, etc.)
      std::string Name_;
   public:
      ///Returns the name of the property
      std::string Name()const{return Name_;}
      ///Returns the MBE truncation order of this property
      int N()const{return Base_t::size();}
      ///Returns the properties for the set of n-mers, indexed by serial number
      PsiMap<LibMolecule::SerialNumber,T>& operator[](const int N){
         return Base_t::operator[](N);
      }
      ///Const version of operator[]
      const PsiMap<LibMolecule::SerialNumber,T>& operator[](const int N)const{
         return Base_t::operator[](N);
      }
      ///Allows direct retrieval of the i-th N-mer's property
      const T& operator()(const int N,const LibMolecule::SerialNumber& i)
         const{return (*this)[N][i];}
      ///Makes a MBEProp object for up to N-mers with that property called Name
      MBEProp<T>(const int N,const std::string& Name="Properties"):
            Base_t(N),Name_(Name){}
      ///Calls += to combine c into the i-th N-mers property
      void Change(const int N,const LibMolecule::SerialNumber& i,const T& c);
      ///Convienent typedef of an iterator of the set of N-mer properties
      typedef typename BaseArg_t::iterator iterator;
      ///const typedef of iterator over set of N-Mer properties
      typedef typename BaseArg_t::const_iterator const_iterator;
      ///Returns an iterator to the first N-mer property
      iterator begin(const int N){return (*this)[N].begin();}
      ///Returns an iterator to one past the last N-Mer property
      iterator end(const int N){return (*this)[N].end();}
      ///Returns a const iterator to the beginning of an N-Mer property
      const_iterator begin(const int N)const{return (*this)[N].begin();}
      ///Returns the const end-point for an N-Mer property
      const_iterator end(const int N)const{return (*this)[N].end();}
      ///Returns a string suited for pretty printing
      std::string PrintOut()const;
};

/* Implemenations below */

template<typename T>
void MBEProp<T>::Change(const int N,const LibMolecule::SerialNumber& i,const T& c){
   if((*this)[N].count(i)==0)(*this)[N][i]=c;
   else (*this)[N][i]+=c;
}

template<typename T>
std::string MBEProp<T>::PrintOut()const{
   std::stringstream Result;
   typedef typename MBEProp<T>::const_iterator It_t;
   for(int i=0;i<N();i++){
      Result<<i+1<<"-body "<<Name_<<":"<<std::endl;
      It_t NMerI=begin(i),NMerEnd=end(i);
      for(;NMerI!=NMerEnd;++NMerI)
         Result<<std::setw(16)<<std::setiosflags(std::ios::fixed)
               <<std::setprecision(16)<<NMerI->second<<std::endl;
   }
   return Result.str();
}

}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_MBEPROP_H_ */
