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
#ifndef SRC_LIB_LIBPARALLEL2_PARALLELENVIRONMENT_H_
#define SRC_LIB_LIBPARALLEL2_PARALLELENVIRONMENT_H_
#include <boost/shared_ptr.hpp>
#include "CommEnvGuts/ParallelEnvironmentGuts.h"
namespace psi{
namespace LibParallel{
class Communicator;

/** \brief This class serves as one of two crucial interfaces to
 *         the new libparallel
 *
 *  The parallel environment is where all the MPI stuff lives.  It is
 *  in charge of keeping track of the current communicator, and giving
 *  communicators out when someone asks.  It also is able to function
 *  if no MPI is occurring.
 *
 *  A program may only contain one parallel environment, hence this class
 *  cannot be copied, or assigned (base class enforces this by having
 *  private copy and assignment operators, as well as no default
 *  constructor).  This class will start MPI, if it hasn't been started,
 *  and it will shut-down MPI when this class is destroyed.  Thus its
 *  destruction needs to be the last thing that occurs.
 *
 */

class ParallelEnvironment{
   private:
      ParallelEnvironmentGuts Guts_;
   public:
      ///Prints Useful debugging info
      void PrintOut()const;

      /** \brief Creates our parallel environment
       *
       *  This should be called the 1st thing made in main.
       *  Older versions of MPI need main's arguments, so we request
       *  them for compatability.
       *
       *  \params[in] argc The number of arguments (including the program
       *                   name), that the program was invoked with
       *  \params[in] argv The list of arguments that the program was
       *                   invoked with.  Length is given by argc.
       */
      ParallelEnvironment(int argc, char* argv[]);

      ///Returns a pointer to the current communicator (you can't change it)
      boost::shared_ptr<const Communicator> GetComm()const;

      ///Expert!!!Returns the original rank of the current MPI process
      int Original()const{return Guts_.Original();}

};
}}//End namespaces



#endif /* SRC_LIB_LIBPARALLEL2_PARALLELENVIRONMENT_H_ */
