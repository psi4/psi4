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

#ifndef SET_H_
#define SET_H_

#include<vector>
#include<iostream>
#include<map>
#include<algorithm>
#include<boost/shared_ptr.hpp>
#include "psi4-dec.h"
namespace psi{
namespace LibFrag{


typedef std::vector<int> SVec;
typedef SVec::iterator VecIt;
typedef VecIt (*Operator)(VecIt, VecIt, VecIt, VecIt, VecIt);

inline static void SetCombine(SVec& ThisSet, const SVec& OtherSet, Operator opr) {
   SVec temp=ThisSet,temp2=OtherSet;
   int MaxSize=ThisSet.size()+OtherSet.size();
   ThisSet.resize(MaxSize); //Make sure there is enough room in ThisSet
   VecIt it=opr(temp.begin(), temp.end(), temp2.begin(), temp2.end(),
         ThisSet.begin());
   ThisSet.resize(it-ThisSet.begin()); //Shrink ThisSet down
}

/** \brief This class is intended to be a mathematical set.
 *
 *  This set works by mapping objects of type T to integers.  You are
 *  free to assign these integers however you like because an std::map is
 *  used for the association.  The objects themselves are not used for the
 *  comparisons and operations, the integers are.  Technically this means
 *  that this set is an index set, and we are performing operations on the
 *  indices.
 *
 *  Because two index sets may correspond to indexing two very different
 *  things care needs to be taken by the programmer that such comparisons
 *  actually make sense.  You are free to make such comparisons as the
 *  set operations are perfectly well defined on such objects.
 *
 */
template <typename T>
class Set{
	private:
        ///Sets are friends regardless of what they contain
        template<typename T2> friend class Set;

		///Function that actually computes the union
        template<typename T2>
        void Union(const Set<T2>& other){
           if(other.size()>1)SetCombine(this->Elements_, other.Elements_,
                 &std::set_union);
        }

		///Function computes the intersection of the two sets
		template<typename T2>
        void Intersection(const Set<T2>& other){
		   if(other.size()>1)SetCombine(this->Elements_, other.Elements_,
		         &std::set_intersection);
		   else this->Clear();  //Intersection with empty set is empty set
		}

		///Fucntion computes the set-difference of this and other
		template<typename T2>
		void SetDiff(const Set<T2>& other){
		   if(other.size()>1)SetCombine(this->Elements_,other.Elements_,
		         &std::set_difference);
		}

		/**\defgroup Relational Implementations
		 *
		 * Note for the following functions one being false does not
         * imply the others are true. If A is not a proper subset of B it is
         * possible that:
         *   A is both a subset and superset of B (they are equal)
         *   A is a proper superset of B, or
         *   A and B have no relation (there are elements of B, not in A).
		 * @{*/
		///Returns true if this is a proper superset of other
		template<typename T2>
		bool IsSuperset(const Set<T2>& other)const{
		   bool is=false;
		   if(other.size()<this->size()){
		      Set<T> temp=other/(*this);
		      if(temp.size()==other.size())is=true;
		   }
		   return is;
		}

		///Returns true if this is a proper subset of other
		template<typename T2>
		bool IsSubset(const Set<T2>& other)const{
		   bool is=false;
		   if(other.IsSuperset(*this))is=true;
		   return is;
		}

		///Returns true if this equals other
		template<typename T2>
		bool IsEqual(const Set<T2>& other)const{
		   bool is=false;
		   if(other.size()==this->size()){
		      Set<T> temp=other/(*this);
		      if(temp.size()==other.size())is=true;
		   }
		   return is;
		}
		/**@}*/

		///Makes this set a copy of other
		void Copy(const Set<T>&other){
		   this->Elements_=other.Elements_;
		   this->Universe_=other.Universe_;
		}

	protected:
	   ///Given an element's number returns that actual element
	   boost::shared_ptr<std::map<int, T> > Universe_;

	/* \brief List of elements in this set
		 *
		 * Elements in this array are assumed to be ordered, from
		 * smallest to largest as determined by int's < operator.
		 * If you add elements to this array
		 * simply call this->Sort() to have the set organize them for
		 * you.
		 *
		 * Note: the C++ defines an object called set, which may
		 * seem appropriate here; however, sets are useful primarily
		 * when you will be searching them and inserting new elements
		 * constantly.  Presently, we will make our sets and then do
		 * operations, likely never inserting again.  Sorted vectors
		 * are better for this purpose.
		 */
		std::vector<int> Elements_;

	public:
		///Returns the actual i-th object
		const T Object(int i)const{return (*Universe_)[i];}

		///Sorts the set in ascending order
		void Sort(){std::sort(Elements_.begin(), Elements_.end());}

		///Default constructor, does nothing
		Set<T>():Universe_(
		      boost::shared_ptr<std::map<int, T> >(new std::map<int,T>())){}

		///Default destructor, no memory to free up
		virtual ~Set<T>(){}

		///Copy constructor.  Guess what it does...
		Set<T>(const Set<T>& other){Copy(other);}

		///Assignment operator, checks for self-assignment
		const Set<T>& operator=(const Set<T>& other){
		   if(this!=&other)Copy(other);return *this;
		}

		///Returns true if this is a *proper* subset of other
		template<typename T2>
		bool operator<(const Set<T2>& other)const{return IsSubset(other);}

		///Returns true if this is a *proper* superset of other
		template<typename T2>
		bool operator>(const Set<T2>& other)const{return IsSuperset(other);}

		///Returns true if this set equals other set
		template<typename T2>
		bool operator==(const Set<T2>& other)const{return IsEqual(other);}

		///Returns true if this set does not equal other set
		template<typename T2>
		bool operator!=(const Set<T2>& other)const{return !((*this)==other);}

		///Returns true if this is a subset (includes equality)
		template<typename T2>
		bool operator<=(const Set<T2>& other)const{
		   return ((*this)==other||(*this)<other);
		}

        ///Returns true if this is a superset (includes equality)
		template<typename T2>
		bool operator>=(const Set<T2>& other)const{
           return ((*this)==other||(*this)>other);
        }

		///Makes this set the union of this set and other set
		template<typename T2>
		void operator*=(const Set<T2>&other){this->Union(other);}

		///Returns the union of this set and another set
		template<typename T2>
		Set<T> operator*(const Set<T2>& other)const{
		   Set<T> Temp(*this);Temp*=other;return Temp;
		}

		///Makes this set the intersection of this set and other set
        template<typename T2>
		void operator/=(const Set<T2>& other){this->Intersection(other);}

		///Returns the intersection of this set and another set
		template<typename T2>
        Set<T> operator/(const Set<T2>& other)const{
		   Set<T> Temp(*this);Temp/=other;return Temp;
		}

		///Makes this set equal to the set difference of this and other
		template<typename T2>
		void operator-=(const Set<T2>& other){
		   this->SetDiff(other);
		}

		///Returns the set-difference of this set and another set
		template<typename T2>
		Set<T> operator-(const Set<T2>& other)const{
		   Set<T> Temp(*this);Temp-=other;return Temp;
		}

		///Adds an atom into this set
		virtual void operator+=(const int NewAtom){
		   this->Elements_.push_back(NewAtom);
		}

        ///Another way of adding atoms
        virtual Set<T> operator<<(const int NewAtom){
           (*this)+=NewAtom;
            return *this;
        }

		///Returns the i-th atom, non-modifiable (0<=i<NAtoms)
		const int& operator[](const int i)const{return Elements_[i];}

		///Returns the i-th atom, modifiable (0<=i<NAtoms)
		//int operator[](const int i){return Elements_[i];}

		///Returns the complement of the current set
		Set<T> Complement()const;

		///Returns the cardinality (size) of this set (lowercase by habit)
		int size()const{return this->Elements_.size();}

		///Prints this set out
		void print_out()const;

		///Check to see if the set contains element i
		bool Contains(const int i)const;

		///Clears the set
		void Clear(){Elements_.clear();}

		///Adds a set to Universe
		virtual void AddSet(int identity,const T& other);


};

template <typename T>
bool Set<T>::Contains(const int i)const{
   SVec temp(1,i);
   SetCombine(temp,this->Elements_,&std::set_intersection);
   return (temp.size()==1);
}

template <typename T>
void Set<T>::print_out()const{
   for(int i=0;i<size();i++)outfile->Printf("%d ",Elements_[i]);
   outfile->Printf("\n");
}

template <typename T>
Set<T> Set<T>::Complement()const{
   Set<T> temp(*this);
   Set<T> temp2();
   std::map<int,T>::iterator itr=Universe_.begin();
   for(; itr!=Universe_.end();++itr)
      temp2<<(itr->first);
   temp2-=temp;
   return temp2;
}

template<typename T>
void Set<T>::AddSet(int identity,const T& other){
   (*Universe_)[identity]=other;
}

}}//End namespaces



#endif /* SET_H_ */
