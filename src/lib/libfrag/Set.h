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
namespace LibFrag{

class Set{
	private:
		///Function that actually computes the union
		void Union(const Set& other);
		///Function computes the intersection of the two sets
		void Intersection(const Set& other);
		///Makes this set a copy of other

		void Copy(const Set&other){
		   this->Atoms=other.Atoms;
		   this->Ghosts=other.Ghosts;
		}

	protected:
		/* \brief List of atoms in this set
		 *
		 * For each atom in this set, j=Atoms[i], means the i-th atom in this set is the j-th atom given in the input.
		 * assumed to be in order (i.e. Atoms[i]<Atoms[j] \forall i<j)*/
		std::vector<int> Atoms;
	public:

        /**For BSSE purposes ghosts attached to this set, made it public for
		calls like Set.Ghosts[i] to distinguish from Set[i], which gives
		the real atom i*/
        std::vector<int> Ghosts;
		void AddGhost(const int i){Ghosts.push_back(i);}

		///Default constructor
		Set(){}
		///Copy constructor
		Set(const Set& other){Copy(other);}
		///Assignment operator
		const Set& operator=(const Set& other){if(this!=&other)Copy(other);return *this;}
		///Makes this set union of this set and other set
		void operator*=(const Set&other){this->Union(other);}
		///Returns the union of this set and another set
		Set operator*(const Set& other){Set Temp(*this);Temp*=other;return Temp;}
		///Makes this set the intersection of this set and other set
		void operator/=(const Set& other){this->Intersection(other);}
		///Returns the intersection of this set and another set
		Set operator/(const Set& other){Set Temp(*this);Temp/=other;return Temp;}
		///Adds an atom into this set
		void operator+=(const int NewAtom){this->Atoms.push_back(NewAtom);}
		///Returns the i-th atom (0<=i<NAtoms)
		const int& operator[](const int i)const{return Atoms[i];}
		int operator[](const int i){return Atoms[i];}
		///Another way of adding atoms
		Set& operator<<(const int NewAtom){(*this)+=NewAtom;return *this;}
		int size()const{return this->Atoms.size();}
		///Prints this set out
		void print(){for(int i=0;i<size();i++)std::cout<<Atoms[i]<<" ";std::cout<<std::endl;}
};
}



#endif /* SET_H_ */
