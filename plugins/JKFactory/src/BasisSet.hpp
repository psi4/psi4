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
#ifndef BASISSET_HPP_
#define BASISSET_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>
#include "JKFactory.hpp"
namespace JKFactory{

class Molecule;
typedef boost::shared_ptr<Molecule> pMol;

/** \brief A Gaussian primitive
 *
 * Gaussian primitives are assumed of the form:
 *
 * \f[
 *    g_i\left(\bf{r}:c_i,\beta_i\right)=c_i\ \rm{e}^{-\beta_i\bf{r}^2},
 * \f]
 *
 * where \f$c_i\f$ is the expansion coefficient associated with this primitive,
 * and \f$\beta_i\f$ is exponent scale factor.  The smaller
 * \f$\beta_i\f$ the more diffuse the primitive.
 */
class Primitive:public JKFactoryBase{
	private:
		///The expansion coefficient of this primitive
		double c;
		///The exponent scale factor of this primitive
		double beta;
	public:
		/**Creates a Primitive with exponent scale factor beta0 and
		expansion coefficient c0*/
		Primitive(double c0=0,double beta0=0):c(c0),beta(beta0){}
		///Returns c
		double GetC()const{return c;}
		///Returns beta
		double GetBeta()const{return beta;}
		///Prints out the primitive
		void PrintOut()const;

};

/** \brief A shell of Gaussian primitives.
 *
 *  Each shell, of momentum \f$\ell\f$, is a set of functions
 *  \f$\left\lbrace\phi^{abc}:a+b+c=\ell\right\rbrace\f$,
 *  each of which is assumed of the form:
 *
 * \f[
 * 		\phi^{abc}\left(\bf{r}\right)=
 * 		\left(x-x_0\right)^a\left(y-y_0\right)^b\left(z-z_0\right)^c
 * 		\sum_i g_i\left(\bf{r}:c_i,\beta_i\right)
 * \f]
 *
 * Note, that the actual member flag controlling whether the shell Cartesian
 * or spherical is true for spherical and false for Cartesian,
 * which is the opposite of how everything outside of the class sees it.  For
 * better or for worse, I decided to do it this way so that I
 * could use words with Cart in them and not interfere with the flag.  Of course I know also have Cartesian coordinates in the class, so
 * that logic is out the window...
 *
 * This class is an abstract base class that is the parent to CartShell and SphereShell, which in turn implement the GetNBasis fxn
 * as is appropriate for that type.
 *
 */
class Shell:public JKFactoryBase{
	private:
		///Returns the letter corresponding to angular moemntum l
		std::string l2string(const int l)const;
	protected:
		///Convenient typedef
		typedef boost::shared_ptr<Primitive> pPrim;
		///True if spherical, false if Cartesian
		bool Spherical;
		///The angular momentum of this shell
		int l;
		///A vector of the Primitives for this shell
		std::vector<pPrim> Gi;
	public:
		/** \brief Creates a Shell
		 *

		 *  \param[in] l_ with angular momentum (default is s-type)
		 *  \param[in] isCart if true, which is the default, orbitals are Cartesian like, if false, they are spherical like
		 */
		Shell(int l_=0,bool IsCart=true):l(l_),Spherical(!IsCart){}
		///Frees memory associated with the shell
		virtual ~Shell(){}
		///Returns the number of primitives
		int GetNPrims()const{return Gi.size();}
		///Returns the angular momentum
		int GetL()const{return l;}
		///Returns true if this a Cartesian like shell
		bool IsCart()const{return (!Spherical);}
		///Returns the number of basis functions in the shell
		virtual int GetNBasis()const=0;
		/** Creates a primitive with expansion coefficient c0 and
		 * exponent scale factor beta0*/
		void AddPrimitive(double c0,double beta0);
		///Returns a pointer to the i-th primitive
		pPrim operator[](const int i){return Gi[i];}
		const pPrim& operator[](const int i)const{return Gi[i];}
		///Prints out the shell
		void PrintOut()const;
};

typedef boost::shared_ptr<Shell> pShell;

/** \brief Specialization of Shell to Cartesian Gaussians
 *
 * In a Cartesian Gaussian of angular momentum \f$\ell=a+b+c\f$, we can assign quantum numbers between 0 and \f$\ell\f$ to each
 * of the three numbers; however, the value of the third is set once we have
 * choosen the first two.  This means we can only really
 * choose two numbers from the possible set.  If we let \f$a=\ell\f$ then
 * \f$b=c=0\f$ necessarily.  If we let \f$a=\ell-1\f$ then
 * \f$\lbrace b=1,c=0\rbrace\f$ or \f$\lbrace b=0,c=1\rbrace\f$.  It is now
 * apparent that in general if \f$a=\ell-\theta\f$, then
 * we have \f$\theta+1\f$ choices of pairs of \f$b\f$ and \f$c\f$, which
 * correspond to the choices for b, which are each value
 * between 0 and \f$\theta\f$.  This means:
 *
 * \f{eqnarray*}{
 *    \rm{NBasis}=&\sum_{i=0}^\ell\left(i+1\right)\\
 *    =&(\ell+1)+\sum_{i=1}^\ell i\\
 *    =&(\ell+1)+\frac{\ell\left(\ell+1\right)}{2}\\
 *    =&\frac{2\ell+2+\ell^2+\ell}{2}\\
 *	  =&\frac{\left(\ell +1\right)\left(\ell +2\right)}{2}
 * \f}
 */
class CartShell:public Shell{
	public:
		///Calls the Shell constructor with IsCart=true
		CartShell(int l0=0):Shell(l0,true){}
		///Returns the number of basis functions
		int GetNBasis()const{return ((GetL()+1)*(GetL()+2)/2);}

};

/** \brief Specialization of Shell to Spherical Gaussians
 *
 * A shell of spherical Gaussians is characterized by a momentum \f$\ell\f$
 * and has one Gaussian for each value in the set
 * \f$\lbrace i:0\le |i|\le\ell\rbrace\f$, which amounts to \f$2\ell+1\f$
 * allowed values.
 *
 */
class SphereShell:public Shell{
	public:
		///Calls Shell constructor with IsCart=false
		SphereShell(int l=0):Shell(l,false){}
		///Returns the number of basis functions
		int GetNBasis()const{return (2*GetL()+1);}
};

/** \brief An abstract base class for Gaussian basis sets
 *
 * The basis set for our system can be thought of as a collection of the basis
 * sets for each atom, each of which is a basis set in its own right.  This class
 * is intended to provide the common functionality for these two basis sets.
 */
class BasisSet:public JKFactoryBase{
	protected:
		///The number of basis functions in this basis
		int NBasis;
		///The number of shells in this basis set
		int NShells;
		///The number of primitives in this basis set
		int NPrims;
		///The maximum momentum
		int max_momentum;
		///The maximum number of primitives
		int max_prims;
		///The ID of the shell with the maximum number of primitives,counting from 0
		int max_prims_ID;
		///The maximum number of basis functions in this basis
		int max_basis;
	public:
		///Sets everything to 0
		BasisSet();
		///No memory allocated=nothing done
		virtual ~BasisSet(){}
		///Returns the number of shells in this basis set
		virtual int GetNShells()=0;
		///Returns the number of basis functions
		virtual int GetNBasis()=0;
		///Returns the maximum momentum
		virtual int GetMaxL()=0;
		///Returns the maximum number of primitives and the ID of the shell with that number
		virtual int GetMaxPrim(int &ID)=0;
		///Returns the maximum number of basis functions in this basis set;
		virtual int GetMaxBasis()=0;
		///Prints the basis set
		virtual void PrintOut()const;
		///Returns the number of primitives
		virtual int GetNPrims();
};

/**** \brief The entire basis set for an atom
 *
 * Not to be confused with the atomic orbital basis set, which is the entire basis set
 * for the molecule, this is the basis set for one atom.
 */
class AtomicBasisSet:public BasisSet{
	private:
		///The coordinates (bohr) of this atom, and therefore the coordinates of the center of each of its shells
		double carts[3];
		///The shells of this atom's basis set
		std::vector<pShell> Shells;
	public:
		/** \brief Creates a basis set for an atom
		 *
		 *  \param[in] Carts The coordinates of the point in space upon which this shell is centered
		 */
		AtomicBasisSet(const double *Carts);
		~AtomicBasisSet(){}
		/** \brief Creates a new shell
		 *
		 * This function creates new shells.  It also keeps track of the angular momentum and the number of shells
		 * \param[in] l The angular momentum of the shell
		 * \param[in] isCart True if the orbitals in this shell are Cartesian-like, for spherical-like use false
		 */
		void AddShell(int l,bool isCart);
		///Returns the number of shells
		int GetNShells(){return NShells;}
		///Returns the number of basis functions
		int GetNBasis();
		///Returns the maximum angular momentum
		int GetMaxL(){return max_momentum;}
		///Returns the maximum number of primitives
		int GetMaxPrim(int &ID);
		///Returns the maximum number of basis functions in a shell
		int GetMaxBasis();
		///Returns the i-th shell of this basis set
		pShell operator[](const int i){return Shells[i];}
		const pShell operator[](const int i)const{return Shells[i];}
		///Prints the basis set
		void PrintOut()const;
};

typedef boost::shared_ptr<AtomicBasisSet> pABS;

/** \brief The top-level basis set object for the JKBuilder library
 *
 * ERD and OED integrals are run by shell, so we need details per shell; however
 * each of those shells retains some properties of the atoms they reside on.  I have
 * made the decision to store the shells for each atom together as objects.  The set of all
 * shells for an atom is a basis set for that atom, an JKBuilder::AtomicBasisSet object.  The
 * set of these basis sets for each atom in a molecular system is then the AOBasisSet.  That is
 * what is stored in this object.
 *
 */
class AOBasisSet:public BasisSet{
	private:
		/** \brief The AO basis
		 *
		 * A vector of the basis sets for each atom, this is the AO basis.  The
		 * order of the AtomicBasisSets that are on it is expected to correspond
		 * to the order the atoms are in the input file.  That is AOBasis[i] should
		 * be the basis set of atom i,counting from 0.
		 */
		std::vector<pABS> AOBasis;
		///Vector containing the number of shells on each atom
		std::vector<int> ShellsPerAtom;
		///Vector containing the number of basis functions per shell
		std::vector<int> BFsPerShell;
	public:
		/** \brief Sets up a basis set for the given molecule
		 *
		 * \param[in] molecule The system we are setting up a basis set for
		 */
		AOBasisSet(pMol molecule);
		///Frees memory
		~AOBasisSet(){}
		///Returns the number of shells in the AObasis
		int GetNShells();
		///Returns the total number of basis functions in the AObasis
		inline int GetNBasis();
		///Returns the maximum momentum contained in any JKBuilder::AtomicBasisSet
		int GetMaxL();
		///Returns the number of atoms this basis is for
		int GetNAtoms(){return AOBasis.size();}
		///Returns the maximum number of primitives in a shell and sets ID to the ID of that shell
		int GetMaxPrim(int &ID);
		///Returns the maximum number of basis functions in a shell
		int GetMaxBasis();
		///Returns a pointer to the i-th atom's basis set
		pABS operator[](const int i){return AOBasis[i];}
		const pABS operator[](const int i)const{return AOBasis[i];}
		///Prints the AO basis set
		void PrintOut()const;
		///Converts shell y, on atom x to an absolute number
		int AbsShellNum(int x,int y);
		///Converts basis fxn z, in shell y, on atom x to an absolute number by means of AbsShellNum
		int AbsBasisNum(int x,int y, int z);
};
}//End namespace JKFactory

#endif /* BASISSET_HPP_ */
