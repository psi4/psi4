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
#ifndef MPITASK_H_
#define MPITASK_H_
namespace psi {
/** \brief Very basic class that stores the properties of an MPI task
 *
 *  For our purposes you need to come up with some labeling scheme for
 *  your tasks.  If it is easier for you to label your tasks with integers
 *  than this object would be of type MPITask<int>, if you want to use
 *  strings, than it is of type MPITask<string>, etc.  i.e.:
 *
 *  \param[in] T The type of your label
 *
 *  Optionally, you many also specify the priority of your tasks.
 *  By default all tasks are assumed to possess the same priority.
 *  Higher numbers mean it needs to be scheduled first, i.e. it is a more
 *  complex task, with complexity being measured relative to other tasks
 *  in your task pool.  It's also worth noting, that currently I only assume
 *  that a task of priority 10 is more complicated than a task of priority 1,
 *  not that it is 10x harder.
 *
 *  If you need to, you may use negative priorities
 *  for tasks that are of an unknown complexity.  For negative priorities
 *  the magnitude is irrelevant, i.e. I do not assume a task of priority -1
 *  is more or less complex than one of priority -2, I simply assume they
 *  both may have a random complexity.  Based on the overall complexity of
 *  your set of tasks an appropriate scheduler will be chosen.
 */
template <typename T>
class MPITask{
   private:
      ///The label the user gave us
      T Label_;

      ///The Priority
      int Priority_;

      ///Deep Copy
      void Copy(const MPITask<T>& other) {
         this->Label_=other.Label_;
         this->Priority_=other.Priority_;
      }

   public:

      ///Makes a task called "Label" that is of complexity "Priority"
      MPITask<T>(const T& Label, const int Priority=0) :
            Priority_(Priority), Label_(Label){}

      ///Makes a new MPITask by copying the old one
      MPITask<T>(const MPITask<T>& other){this->Copy(other);}

      ///Assignment operator
      const MPITask<T>& operator=(const MPITask<T>& other){
         if (this!=&other) this->Copy(other);
         return *this;
      }

      ///Returns the label
      T GetLabel() const {return Label_;}

      ///Returns the priority
      int Priority()const{return Priority_;}

      ///No memory to free up
      ~MPITask<T>() {}
};
}//End namespace psi
#endif /* MPITASK_H_ */
