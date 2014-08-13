/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include "JKFactory.hpp"
#include <vector>
namespace JKFactory{
class Atom:public JKFactoryBase{
	private:
		///Atomic number (between 1 and 36 for now)
		int Z;
		///The number of electrons
		int NElectrons;
		///Carts in bohr
		double Carts[3];
		///Function that returns the chemical symbol of atomic number Z in caps
		std::string ZtoAt(const int Z) const;
	public:
		/** \brief Creates an atom
		 *
		 * \param[in] Z_ is the atomic number
		 * \param[in] carts_ is a 3 element array where carts_[0],carts[1], and carts[2] are the x,y, and z coordinates of the atom
		 * \param[in] bohr True means the carts that are passed in are in atomic units, false means they are in Angstroms, default is false
		 */
		Atom(const int Z_,const double *Carts_);
		~Atom(){}
		///Returns the i-th coordinate of the atom (x=0,y=1,z=2)
		const double& operator[](const int i)const{return Carts[i];}
		///Returns the atomic number
		int GetZ()const{return Z;}
		///Prints the Atom class out for debugging purposes
		void PrintOut()const;
		///Returns the number of electrons (calls GetZ)
		int NElec()const{return NElectrons;}
};

class Molecule:public JKFactoryBase{
	private:
		///Convent typedef of
		typedef boost::shared_ptr<Atom> pAtom;
		///The atoms in this molecule
		std::vector<pAtom> Atoms;
		///The total number of electrons (adjusted for charge)
		int NElectrons;
	public:
		///Lists the atoms in the molecule
		void PrintOut()const;
		int NElec()const{return NElectrons;}
		///Returns the number of atoms
		int NAtoms()const{return Atoms.size();}
		///Short-cut to the coordinates of the atoms
		const double& operator()(const int atom,const int comp)const{
			return (*(*this)[atom])[comp];
		}
		///Returns a shared pointer to atom i
		const pAtom operator[](const int i) const{return Atoms[i];}
		/** Given an Natoms long vector Zs of atomic numbers and a 3*Natoms
		 *  long vector of Cartesian coordinates (in bohr) makes a new
		 *  Molecule class
		 */
		Molecule(const std::vector<int>& Zs,const std::vector<double>& Carts,
				const int charge=0);
		///No memory allocated (under our control at least)
		~Molecule(){}

};

}
#endif /* MOLECULE_HPP_ */
