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
#ifndef TOOLSET_H_
#define TOOLSET_H_

#include <boost/shared_ptr.hpp>
#include <vector>
#include "../LibParallelBase.h"
#include "../TaskJobGuts/MPITaskGuts.h"

namespace psi{
namespace LibParallel{
class Communicator;
/** \brief This class is intended to be the base class for operations between
 *         Masters and Slaves
 *
 *
 */
class ToolSet: public LibParallelBase{
   protected:
      enum MessageCodes{ANY=-1,NEXT=1,TASKSIZE=2,TASKS=3};
      boost::shared_ptr<const Communicator> State_;
   public:
      ToolSet(boost::shared_ptr<const Communicator>& State);
      virtual ~ToolSet();
      virtual bool Done()=0;
      virtual void GiveTasks()=0;
      virtual void PrintOut()const=0;
};


}}



#endif /* TOOLSET_H_ */
