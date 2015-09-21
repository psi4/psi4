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
#ifndef PSIAPI_UTILS_PSISET_H_
#define PSIAPI_UTILS_PSISET_H_

namespace PsiAPI{

/** \brief A generic set class
 *
 *   In my coding I quite often want a true set, and possibly
 *   built out of a container other than std::set.  This class
 *   inherits from a container of your choice, giving it the
 *   same basic functions.  Then on top of those functions
 *   I define union, intersection, set-difference, and symmetric-difference.
 *
 *   To use this class your base container must be weakly sorted and
 *   the elements of your container must have a less than operator
 *   defined.  I'm lazy and use the std::algorithms, which do not
 *   occur in place, thus all operators incur a copy
 *
 */
template<typename T>
class PsiSet: public T{
   public:
      ///Union
      const PsiSet<T>& operator*=(const PsiSet<T>& other);
      PsiSet<T> operator*(const PsiSet<T>& other)const;
      ///Intersection
      const PsiSet<T>& operator/=(const PsiSet<T>& other);
      PsiSet<T> operator/(const PsiSet<T>& other)const;
      ///Symmetric difference
      const PsiSet<T>& operator+=(const PsiSet<T>& other);
      PsiSet<T> operator+(const PsiSet<T>& other)const;
      ///Set difference
      const PsiSet<T>& operator-=(const PsiSet<T>& other);
      PsiSet<T> operator-(const PsiSet<T>& other)const;
};

/****************** Implementations *********************/
template<typename T>
const PsiSet<T>& PsiSet<T>::operator*=(const PsiSet<T>& other){
   PsiSet<T> Temp;
   std::set_union(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   (*this)=Temp;
   return *this;
}

template<typename T>
const PsiSet<T>& PsiSet<T>::operator/=(const PsiSet<T>& other){
   PsiSet<T> Temp;
   std::set_intersection(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   (*this)=Temp;
   return *this;
}

template<typename T>
const PsiSet<T>& PsiSet<T>::operator+=(const PsiSet<T>& other){
   PsiSet<T> Temp;
   std::set_symmetric_difference(this->begin(),this->end(),
         other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   (*this)=Temp;
   return *this;
}

template<typename T>
const PsiSet<T>& PsiSet<T>::operator-=(const PsiSet<T>& other){
   PsiSet<T> Temp;
   std::set_difference(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   (*this)=Temp;
   return *this;
}

template<typename T>
PsiSet<T> PsiSet<T>::operator*(const PsiSet<T>& other)const{
   PsiSet<T> Temp;
   std::set_union(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   return Temp;
}

template<typename T>
PsiSet<T> PsiSet<T>::operator/(const PsiSet<T>& other)const{
   PsiSet<T> Temp;
   std::set_intersection(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   return Temp;
}

template<typename T>
PsiSet<T> PsiSet<T>::operator+(const PsiSet<T>& other)const{
   PsiSet<T> Temp;
   std::set_symmetric_difference(this->begin(),this->end(),
         other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   return Temp;
}

template<typename T>
PsiSet<T> PsiSet<T>::operator-(const PsiSet<T>& other)const{
   PsiSet<T> Temp;
   std::set_difference(this->begin(),this->end(),other.begin(),other.end(),
         std::inserter(Temp,Temp.begin()));
   return Temp;
}

}
#endif /* PSIAPI_UTILS_PSISET_H_ */
