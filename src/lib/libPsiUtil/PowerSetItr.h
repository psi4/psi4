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
#ifndef SRC_LIB_LIBPSIUTIL_POWERSETITR_H_
#define SRC_LIB_LIBPSIUTIL_POWERSETITR_H_

namespace psi{

/** \brief Class to facilitate iterating over the power set of a set
 *
 *   Again I like generality so this iterator will work with any
 *   container that defines the operations insert, assignment, and
 *   const_iterator, which should be all containers in the std template
 *   library
 *
 *   Assume we have some set \f$S\f$ that has \f$N\f$ elements in it.
 *   The powerset of \f$S\f$, denoted \f$\mathbb{P}(S)\f$ is the set:
 *   \f$\mathbb{P}(S)=\lbrace s_i\ :\ s_i\ \subseteq\ S\rbrace\f$.  Defining
 *   \f$p_0=\emptyset\f$ we iteratively form $\mathbb{P}(S)\f$ in \f$N\f$
 *   steps, where in step \f$i\f$ we form the set:
 *   \f$p_i=p_{i-1}\bigcup_{j\in p_{i-1}}\lbrace  s_i\ \cup\ j\rbrace\f$
 *
 *   Assuming you gave the set \f$S={i, j, k, l}\f$, the results would
 *   be:
 *   \verbatim
 *   (empty set)
 *   i
 *   j
 *   i j
 *   k
 *   i k
 *   j k
 *   i j k
 *   l
 *   i l
 *   j l
 *   i j l
 *   k l
 *   i k l
 *   j k l
 *   i j k l
 *   \endverbatim
 *   Another way of thinking about this is that they follow the order
 *   of counting in binary from 0 to the number of elements in your set
 *   (this is actually how I implemented ).
 *
 *   \param T The type of your set (e.g. std::vector<double>)
 *
 */
template<typename T>
class PowerSetItr{
   private:
      typedef std::vector<bool> Set_t;
      const T Orig_;
      Set_t InSet_;
      T Current_;
      ///Keeps track of how far to the left in InSet we actually got
      size_t Edge_;
      ///Counts in binary using InSet_ (calls Fill if we aren't done)
      void Next();
      ///Sets current off of InSet_
      void Fill();
   public:
      PowerSetItr(const T& Set):
         Orig_(Set),InSet_(Set.size(),false),Edge_(0){}
      const T& operator*()const{return Current_;}
      const T* operator->()const{return &Current_;}
      bool Done()const{return Edge_==InSet_.size()+1;}
      const PowerSetItr& operator++(){Next();return *this;}
};

/***************Implementation****************/

template<typename T>
void PowerSetItr<T>::Next(){
   int size=InSet_.size();
   for(int i=size-1;i>=0;--i){
      //Find first element, closest to the right that is false
      if(!InSet_[i]){
         if(i+Edge_==size-1)Edge_++;
         //Negate it (now a 1), and all the elements (all 1's) right of it
         for(;i<size;++i)InSet_[i]=!InSet_[i];
      }
      if(i==size){
         Fill();
         break;
      }
      //Make sure we trip the iterator termination condition
      else if(i==0)Edge_++;
   }
}

template<typename T>
void PowerSetItr<T>::Fill(){
   Set_t::reverse_iterator It=InSet_.rbegin();
   T Temp;
   typename T::const_iterator ElemI=Orig_.begin();
   for(size_t counter=0;counter<Edge_;++It,++ElemI,++counter)
      if(*It)Temp.insert(Temp.end(),*ElemI);
   Current_=Temp;

}

}//End namespace



#endif /* SRC_LIB_LIBPSIUTIL_POWERSETITR_H_ */
