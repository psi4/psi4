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
#ifndef TASKSTATISTICS_H_
#define TASKSTATISTICS_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include "../Algorithms.h"
#include "../Communicator.h"
#include "../LibParallelBase.h"
#include "../libPsiUtil/PsiMap.h"
namespace psi{
namespace LibParallel{
class MPITaskGuts;
/** \brief Class for determining the properties of a task set
 *
 *
 *  This is really the class that controls which scheduling algorithm
 *  is chosen.  In order to choose the algorithm we first see if there
 *  are any unknown priorities.  If there are, I simply assign a random
 *  priority to each task and flow continues to the next part.  In the
 *  future it may be better to consider a less drastic approach, such
 *  as scheduling around the few unknown tasks.
 *
 *  Now with the priority of each task known we determine how many different
 *  priorities we have, and how many of each.  If we only have one priority
 *  we make a simple scheduler which assigns tasks: process*batchisize through
 *  (process+1)*batchsize to each process.
 *
 *  If we assigned random priorities we always use a dynamic scheduler.
 *
 *
 *  If you have requested a dynamic scheduler, you will
 *  always get one, and as of right now that scheduler will be of
 *  type DynamicRoundRobin.  If you have a roughly even distribution of
 *  each priority you will get a static round robin scheduler
 *
 */
class TaskStats: public LibParallelBase{
   private:
      ///The suggested algorithm
      SchedAlgorithm Sugg_;

      ///The number of times each priority appears
      PsiMap<int,int> PriBins_;

      ///Whether there are negative priority tasks or not
      bool UnknownPri_;

      ///Handles unknown tasks as described in class description
      void HandleUnknown(
            boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks,
            boost::shared_ptr<const Communicator>& State);

   public:
      ///Returns the number of priorities
      int NPriorities()const;

      ///Returns the number of times a given priority occurs
      int NTimes(const int i)const;

      ///The recommended algorithm
      SchedAlgorithm SuggestedAlgorithm()const;

      /** \brief Sets up the TaskStats object
       *
       *
       *   \param[in,out] Tasks The set of tasks the stats are for
       *                        Unknown complexity will be removed
       *   \param[in] ForceDyn Should we force the object to pick
       *                       a dynamic scheduler?
       */
      TaskStats(boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks,
            const bool ForceDyn,boost::shared_ptr<const Communicator>& State);

};


}}//End namespace



#endif /* TASKSTATISTICS_H_ */