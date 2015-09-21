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
#ifndef PSIAPI_UTILS_NESTEDITERATOR_H_
#define PSIAPI_UTILS_NESTEDITERATOR_H_
#include<utility>
namespace PsiAPI{

/** \brief A general iterator class for a container that contains a container
 *
 *  When making hierarchal classes one sometimes wants the ability to
 *  loop over a lower class from a higher one.  For example in a Gaussian
 *  basis set the actual basis set is made up of shells which in turn
 *  are made up of primitives.  If we wanted to loop over primitives we
 *  would first loop over orbitals, then for each orbital we would loop
 *  over it's primitives.  This class wraps that dual iteration into one
 *  iterator
 *
 *  \param[in] T1 The type of the parent iterator
 *  \param[in] T2 The type of the iterator contained inside the parent
 *  \param[in] R  The type of the object returned upon dereferencing this
 *                iterator
 *
 */
template<typename T1, typename T2, typename R>
class NestedIteratorBase{
   protected:
      ///The type of the parent iterator set
      typedef typename std::pair<T1,T1> Parent_t;
      ///The type of the child iterator set
      typedef typename std::pair<T2,T2> Child_t;

      ///The current and end parent iterators
      Parent_t ParentItr_;
      ///The sub iterators of the parent iterator
      Child_t ChildItr_;
   public:
      ///\name Constructor and Destructors
      ///@{
      ///Member functions designed to make NestedIterators

      ///Makes a new iterator (Begin==End) gives termination condition
      NestedIteratorBase(){}

      ///No memory to free up
      virtual ~NestedIteratorBase(){}
      ///@}

      ///Forward prefix-increment the nested iterator
      const NestedIteratorBase<T1,T2,R>& operator++();

      ///Returns true if two iterators are equal
      bool operator==(const NestedIteratorBase<T1,T2,R>& other)const;

      ///Returns true if two operators are not equal
      bool operator!=(const NestedIteratorBase<T1,T2,R>& other)const{
         return !(operator==(other));
      }

      virtual R& operator*()=0;
      virtual const R& operator*()const=0;
      R* operator->(){return &operator*();}
      const R* operator->()const{return &operator*();}
};

///Specialization to normal iterators
template<typename T1,typename T2, typename R>
class NestedIterator:public NestedIteratorBase<T1,T2,R>{
   public:
      NestedIterator(T1 BeginItr,T1 EndItr):ParentItr_(BeginItr,EndItr){
         if(BeginItr!=EndItr)
            this->ChildItr_=
                 Child_t(BeginItr->begin(),BeginItr->end());
      }
      R& operator*(){return *(this->ChildItr_.first);}
      const R& operator*()const{return *(this->ChildItr_.first);}
};

///Slight specialization to when the nested array contains pointers
template<typename T1,typename T2,typename R>
class NestedPtrItr: public NestedIteratorBase<T1,T2,R>{
   public:
      NestedPtrItr(T1 BeginItr,T1 EndItr){
         if(BeginItr!=EndItr)
           this->ChildItr_=
                 Child_t((*BeginItr)->begin(),(*BeginItr)->end());
      }
      R& operator*(){return *(*(this->ChildItr_.first));}
      const R& operator*()const{return *(*(this->ChildItr_.first));}
};

/************************ Implementations ****************/

template<typename T1,typename T2,typename R>
const NestedIteratorBase<T1,T2,R>& NestedIteratorBase<T1,T2,R>::operator++(){
   if(ChildItr_.first!=ChildItr_.second)//Can increment child
      ++(ChildItr_.first);
   else{
      ++(ParentItr_.first);
      if(ParentItr_.first!=ParentItr_.second)
         ChildItr_=
               Child_t(ParentItr_.first->begin(),ParentItr_.first->end());
   }
}
template<typename T1,typename T2,typename R>
bool NestedIteratorBase<T1,T2,R>::
operator==(const NestedIteratorBase<T1,T2,R>& other)const{
   if(ParentItr_.first!=other.ParentItr_.first)return false;
   else if(ParentItr_.first!=ParentItr_.second)
      return (ChildItr_.first==other.ChildItr_.first);
   return true;
}

}//End PsiAPI namespace
#endif /* PSIAPI_UTILS_NESTEDITERATOR_H_ */
