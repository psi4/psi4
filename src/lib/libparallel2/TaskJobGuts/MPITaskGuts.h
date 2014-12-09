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
#ifndef MPITASKGUTS_H_
#define MPITASKGUTS_H_
#include <utility>
#include "../LibParallelBase.h"
namespace psi{
namespace LibParallel{

class MPITaskGuts:public LibParallelBase{
   private:
      ///The actual task this object is wrapping
      std::pair<int,int> Task_;

      ///Makes a deep copy
      void Copy(const MPITaskGuts& other);
   public:
      ///Sets first element of Task_ to Priority, 2nd to Number
      MPITaskGuts(const int Priority,const int Number);

      ///Creates this MPITaskGuts as a copy of other
      MPITaskGuts(const MPITaskGuts& other);

      ///Assigns copy of other to this
      const MPITaskGuts& operator=(const MPITaskGuts& other);

      ///True if this task has a higher priority (its Task_.first is less)
      bool operator<(const MPITaskGuts& other)const;

      ///True if this task is less than or equal to other
      bool operator<=(const MPITaskGuts& other)const;

      ///True if this task is greater than other
      bool operator>(const MPITaskGuts& other)const;


      ///True if this task is greater or equal
      bool operator>=(const MPITaskGuts& other)const;

      ///Equal if two tasks are of the same priority
      bool operator==(const MPITaskGuts& other)const;

      ///Returns the Priority, 0 is highest priority
      int Priority()const;

      ///Allows changes to priority
      void SetPriority(const int i);

      ///Returns the Number
      int Number()const;

      ///Debugging function that prints out the task
      void PrintOut()const;
};


}}//End namespaces



#endif /* MPITASKGUTS_H_ */
