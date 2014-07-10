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
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef JKFACTORY_HPP_
#define JKFACTORY_HPP_

#include<sstream>
#include<boost/shared_ptr.hpp>

namespace JKFactory{
class MPIManager;

///The base class for the JKFactory library
class JKFactoryBase{
	private:
		std::stringstream buffer;
	protected:
		///The object that manages the MPI Environment
		static boost::shared_ptr<MPIManager> MPI;
	public:
		///Prints an error message and then crashes the program
		void Error(const std::string& message)const{Print(message);exit(1);}
		///provides a uniform method for printing throughout the library
		template<class T>
		void operator<<(T& message)const{
		      std::stringstream temp;
		      temp<<message;
		      Print(temp.str());
		}
		void Print(const std::string& message)const{std::cout<<message;}
		///All classes should override this function so that it print's out the class
		virtual void PrintOut()const=0;
		///Destructor does notin' 'cause there ain't anythin' to do
		virtual ~JKFactoryBase(){}
};

}



#endif /* JKFACTORY_HPP_ */
