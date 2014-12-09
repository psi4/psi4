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
#ifndef SRC_LIB_LIBPARALLEL2_COMMENVGUTS_PARALLELENVIRONMENTGUTS_H_
#define SRC_LIB_LIBPARALLEL2_COMMENVGUTS_PARALLELENVIRONMENTGUTS_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "../Communicator.h"
#include "../LibParallelBase.h"

///Forward declaration of boost's environment
namespace boost{
namespace mpi{
class environment;
}}

namespace psi{
namespace LibParallel{
class ParallelEnvironment;
class ParallelEnvironmentGuts: public LibParallelBase{
   private:
      ///No copying
      ParallelEnvironmentGuts(const ParallelEnvironment& other){}

      ///No assignment
      const ParallelEnvironmentGuts& operator=(const ParallelEnvironmentGuts& other)
            {return *this;}
      friend class ParallelEnvironment;
   protected:
      ///The actual environment (only made if MPI is present)
      boost::shared_ptr<boost::mpi::environment> Env_;

      ///Communicators
      std::vector<boost::shared_ptr<Communicator> > Comms_;

      ParallelEnvironmentGuts(int argc, char* argv[]);
   public:
      virtual ~ParallelEnvironmentGuts(){}
      virtual void PrintOut()const;
      void UpdateComms();
      int Original()const;
      void AddComm(boost::shared_ptr<Communicator> Comm);
      virtual boost::shared_ptr<const Communicator> GetComm()const;
};

}}//End namespaces


#endif /* SRC_LIB_LIBPARALLEL2_COMMENVGUTS_PARALLELENVIRONMENTGUTS_H_ */
