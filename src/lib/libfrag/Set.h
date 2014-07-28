/*
 * Set.h
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */

#ifndef SET_H_
#define SET_H_

#include<vector>
#include<iostream>
#include<map>
namespace LibFrag{

/** \brief This class is intended to be a mathematical set, although I call
 *  the elements atoms, there is no chemical level to them in this class.
 */

class Set{
	private:
		///Function that actually computes the union
		void Union(const Set& other);
		///Function computes the intersection of the two sets
		void Intersection(const Set& other);

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
		bool IsSuperset(const Set& other)const;
		///Returns true if this is a proper subset of other
		bool IsSubset(const Set& other)const;
		///Returns true if this equals other
		bool IsEqual(const Set& other)const;
		/**@}*/

		///Makes this set a copy of other
		void Copy(const Set&other){
		   this->Atoms=other.Atoms;
		}

	protected:

	/* \brief List of atoms in this set
		 *
		 * For each atom in this set, j=Atoms[i], means the i-th atom in this set is the j-th atom given in the input.
		 * assumed to be in order (i.e. Atoms[i]<Atoms[j] \forall i<j)*/
		std::vector<int> Atoms;

	public:

		///Sorts the set
		void Sort();

		///Default constructor
		Set(){}

		///Default destructor
		virtual ~Set(){}

		///Copy constructor
		Set(const Set& other){Copy(other);}

		///Assignment operator
		const Set& operator=(const Set& other){
		   if(this!=&other)Copy(other);return *this;
		}

		///Returns true if this is a *proper* subset of other
		bool operator<(const Set& other)const{return IsSubset(other);}

		///Returns true if this is a *proper* superset of other
		bool operator>(const Set& other)const{return IsSuperset(other);}

		///Returns true if this set equals other set
		bool operator==(const Set& other)const{return IsEqual(other);}

		///Returns true if this is a subset (includes equality)
		bool operator<=(const Set& other)const{
		   return ((*this)==other||(*this)<other);
		}

        ///Returns true if this is a superset (includes equality)
        bool operator>=(const Set& other)const{
           return ((*this)==other||(*this)>other);
        }

		///Makes this set union of this set and other set
		void operator*=(const Set&other){this->Union(other);}

		///Returns the union of this set and another set
		Set operator*(const Set& other)const{
		   Set Temp(*this);Temp*=other;return Temp;
		}

		///Makes this set the intersection of this set and other set
		void operator/=(const Set& other){this->Intersection(other);}

		///Returns the intersection of this set and another set
		Set operator/(const Set& other)const{
		   Set Temp(*this);Temp/=other;return Temp;
		}

		///Adds an atom into this set
		void operator+=(const int NewAtom){this->Atoms.push_back(NewAtom);}

		///Returns the i-th atom, non-modifiable (0<=i<NAtoms)
		const int& operator[](const int i)const{return Atoms[i];}

		///Returns the i-th atom, modifiable (0<=i<NAtoms)
		int operator[](const int i){return Atoms[i];}

		///Another way of adding atoms
		Set& operator<<(const int NewAtom){(*this)+=NewAtom;return *this;}

		///Returns the cardinality (size) of this set
		int size()const{return this->Atoms.size();}

		///Prints this set out
		void print_out()const;

		///Check to see if the set contains element i
		bool contains(const int i);
};


}



#endif /* SET_H_ */
