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
#include <iomanip>
#include "exception.h"
namespace psi{
namespace LibFrag{

/** \brief A class to hold generic Many-Body properties
 *
 *   The template parameter T is the type of the property.
 *   If your property defines:
 *   \code
 *   const T& operator+=(const T& other);
 *   const T& operator*(const double scalar)const;
 *   \endcode
 *
 *   then the MBE and GMBE will be able to reassemble your property,
 *   automatically for you.
 *
 *
 */

template<typename T>
class MBEProp:protected std::vector<std::vector<T> >{
   private:
      typedef std::vector<std::vector<T> > Base_t;
      std::string Name_;
   public:
      int N()const{return Base_t::size();}
      std::vector<T>& operator[](const int N){
         return Base_t::operator[](N);
      }
      const std::vector<T>& operator[](const int N)const{
         return Base_t::operator[](N);
      }
      const T& operator()(const int N,const int i)const{return (*this)[N][i];}
      MBEProp<T>(const int N,const std::string& Name="Properties"):
            Base_t(N),Name_(Name){}
      void Change(const int N,const int i,const T& c);
      typedef typename std::vector<T>::iterator iterator;
      typedef typename std::vector<T>::const_iterator const_iterator;
      iterator begin(const int N){return (*this)[N].begin();}
      iterator end(const int N){return (*this)[N].end();}
      const_iterator begin(const int N)const{return (*this)[N].begin();}
      const_iterator end(const int N)const{return (*this)[N].end();}
      std::string PrintOut()const;
};

/* Implemenations below */

template<typename T>
void MBEProp<T>::Change(const int N,const int i,const T& c){
   int Size=(*this)[N].size();
   if(Size<i)throw PSIEXCEPTION("Not Enough Values");
   else if(Size==i)(*this)[N].push_back(c);
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
         Result<<std::setiosflags(std::ios::fixed)
               <<std::setprecision(16)<<(*NMerI)<<std::endl;
   }
   return Result.str();
}




}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_MBEPROP_H_ */
