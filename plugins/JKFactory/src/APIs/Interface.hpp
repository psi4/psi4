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

#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <mpi.h>
#include "GTFockInter.hpp"
#include "MPIManager.hpp"
namespace JKFactory{
class Molecule;
class AOBasisSet;
///Convenient typedef for the Basis Set
typedef boost::shared_ptr<JKFactory::AOBasisSet> pMyBasis;
///Convenient typedef for the Molecule
typedef boost::shared_ptr<JKFactory::Molecule> pMol;
/** \brief The base class from which all APIs should derive
 *
 *  To add support for a new electronic structure package (ESP)
 *  derive a class from this class.  The resulting class will be called "Your
 *  Interface".  Then have your ESP make an instance of "Your Interface", not
 *  this one.  "Your interface" should take data from your ESP and use it to
 *  fill in the Molecule and BasisSet objects appropriately.  Follow the
 *  instructions found in those class descriptions for more details.  It is
 *  important that these two classes are set up first as the remainder of the
 *  set-up will look to these two classes for parameters.
 *
 *  Next you will have to set up the Js and Ks vectors.  These are designed
 *  to be wrappers around the memory location that you will be using to store
 *  the results of the calculations.  This library will not allocate the
 *  memory for those locations, nor will it free that memory.
 *
 *  For example, assume I have an SCF routine managed by a class SCF:
 *
 *  \code
 *  class SCF{
 *     double* J;
 *     double* K;
 *     std::vector<double> ManyJs;
 *     std::vector<double> ManyKs;
 *  };
 *  \endcode
 *
 *  When you create the instance of "Your Interface" it would be something like
 *  this (top for a single J and K, bottom for multiples) :
 *
 *  \code
 *  //One J and K matrix
 *  UpdateJ(SCF::J);
 *  UpdateK(SCF::K);
 *  //More than one J and K matrix
 *  UpdateJ(SCF::ManyJs);
 *  UpdateK(SCF::ManyKs);
 *  \endcode
 *
 *  The number of J matrices needs to equal the number of K matrices, unless
 *  you do not want K matrices, then K must be 0.  That is the set-up phase.
 *
 *  When you want to actually compute J and K call Interface::UpdateDensity(X),
 *  where X is a std::vector<double*> of density matrices if you want multiple
 *  Js and Ks.  If you want only a single J and K call the UpdateDensity where
 *  X is a double*.  After the call, Js and Ks will contain your new J and K
 *  matrices.
 *
 */
class Interface:private GTFockInter{
	protected:


		///The molecular system that the density corresponds to
		pMol System;
		///The basis set for the Molecule contained in System
		pMyBasis Basis;

	public:
		/** Calling this UpdateDensity with a vector of density matrices that
		 * is N matrices long.  Will cause Js and Ks to be populated with N
		 * J matrices and N K matrices respectively, with Js[0] corresponding
		 * to NewDensitys[0], etc.
		 */
		void UpdateDensity(std::vector<double*>& NewDensitys);
		///Same as above UpdateDensity except for a single density matrix
		void UpdateDensity(double* NewDensity);

		///Updates Js so that its doubles point to those in NewJs
		void UpdateJ(std::vector<double*>& NewJs);
		///Wrapper to other UpdateJ, for one J
		void UpdateJ(double* J);

		///Updates Ks, so that its doubles point to those in NewKs
		void UpdateK(std::vector<double*>& NewKs);
		///Wrapper to other UpdateK, for one K
		void UpdateK(double* K);

		/**Makes a new interface, needs to know the integral screening
		 * threshold and optionally whether or not the density matrices are
		 * symmetric.  All matrices are assumed symmetric unless you specify
		 * otherwise.
		 */
		Interface(double ScreenThresh,const MPI_Comm& CurrentComm,
		      bool Symm=true):
			GTFockInter(ScreenThresh,Symm){MPI->Initialize(CurrentComm);
			std::string Comm=MPI->Comm();
			    if(MPI->Me(Comm)==0){
			       std::cout<<"Starting JKFactory!!!!"<<std::endl;
			       std::cout<<"JKFactory is using "<<MPI->NProc(Comm)<<" processes.";
			       std::cout<<std::endl;
			    }
		}

		///Shared pointers handle memory deallocation, so nothing to do
		virtual ~Interface(){}

		///Useful debugging function that will print lots of details to screen
		virtual void PrintOut()const;



};

}



#endif /* INTERFACE_HPP_ */
