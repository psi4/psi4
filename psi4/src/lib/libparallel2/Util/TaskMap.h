/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#ifndef TASKMAP_H_
#define TASKMAP_H_

#include "../../libPsiUtil/PsiMap.h"
#include <vector>
#include <iostream>
#include "../LibParallelBase.h"
namespace psi{
namespace LibParallel{

/** \brief A map to help me keep my sanity
 *
 *  I define the absolute task order as the order the user gave them to
 *  me in.  The relative task order is the order they were actually run
 *  in.  This class is mainly designed to be used in conjunction with
 *  putting a user's data back in order.
 */

class TaskMap: public LibParallelBase{
   private:
      ///Map from relative task number 2 slave ID
      PsiMap<int,int> Task2Slave_;

      ///Map from relative task 2 absolute
      PsiMap<int,int> Rel2Abs_;

      ///Map from absolute 2 relative
      PsiMap<int,int> Abs2Rel_;

      ///The number of tasks each slave has
      std::vector<int> TasksPerSlave_;

      ///TaskQueue_[i][j] is the j-th task slave i ran
      std::vector<std::vector<int> > TaskQueue_;

   public:

      ///Makes a task map capable of monitoring work distributions to NSlaves
      TaskMap(const int NSlaves);

      ///Returns the slave associated with a given Absolute task
      int Slave(const int Abs)const;

      ///Returns the relative index, given an absolute
      int Rel(const int Abs)const;

      ///Returns the absolute index, given a relative one
      int Abs(const int Rel)const;

      ///Returns the index where data regarding Absolute task Abs resides
      int Index(const int Abs,const int slave=-1)const;

      void AddTask(const int Rel,const int Abs,const int slave);

      ///Returns the offset (in tasks) for the start of slave "slave"'s task
      int Offset(const int slave)const;

};


}}//End namespaces



#endif /* TASKMAP_H_ */