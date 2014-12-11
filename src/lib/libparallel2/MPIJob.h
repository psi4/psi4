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
#ifndef MPIJOB_H_
#define MPIJOB_H_

#include "TaskJobGuts/MPIJobGuts.h"
#include "Algorithms.h"
namespace psi {

/** \brief The object that handles all your MPI needs, designed to be the
 *         public interface.  The implementation is contained in the
 *         MPIJobGuts class for the interested parties.
 *
 *  In getting ready to work with this class, it helps if we first define
 *  some key terms.  To facilitate this, let  us work with a concrete
 *  parallel problem, namely let's assume we are parallelizing the
 *  calculation of an interaction energy between a water molecule, and a
 *  chloride ion, by giving each of the 3 calculations to it's own MPI
 *  processes.  Then we define:
 *
 *  - Task:the series of commands that we are actually performing in parallel.
 *    - For our interaction energy we have three tasks, computing the
 *      energy of:
 *        - The Complex,
 *        - The chloride ion, and
 *        - The water molecule
 *  - Priority=The relative complexity of each task, used for determining the
 *    order tasks will be run in, and for load balancing purposes.
 *    - Priority is assigned assuming 0 is the most complex
 *    - Negative priorities mean you have no idea how complex a task is
 *    - For our interaction energy we have three priorities:
 *       - 0 for the complex,
 *       - 1 for the chloride ion (assuming it has more basis fxns
 *         than water), and
 *       - 2 for the water molecule
 *
 *  Really, that's pretty much all you need to know theory-wise, so
 *  switching over to how to use these objects.  When it comes time
 *  to perform your MPI tasks, you first create a vector of
 *  boost::shared_ptrs to the MPITasks you want to perform.  Again, the
 *  label is for your purposes, so make it convenient for you.  This is also
 *  the design philosophy behind using shared_ptr's, namely to avoid copying
 *  a potentially complicated label.
 *
 *  Assuming you have N tasks, all of which are approximately the same
 *  priority, and are easily mappable to integers, a typical usage case
 *  would be something like:
 *  \code
 *  //Typedefs to avoid typing complicated type names all over
 *  typedef MPITask<int> TaskType;
 *  typedef boost::shared_ptr<TaskType> SharedTask;
 *
 *  std::vector<SharedTask> TaskQue;
 *  for(int TaskNumber=0;TaskNumber<N;TaskNumber++){
 *      int priority=0;
 *      SharedTask TempPtr(new TaskType(TaskNumber,priority));
 *      TaskQue.push_back(TempPtr);
 *  }
 *  \endcode
 *
 *  Once you have the array of tasks, you then create an MPI job, which
 *  is creatively, an object with name MPIJob.  You then loop over
 *  the tasks it gives you.  The following code shows you how to do this
 *  in one of two types of loop: a for loop or a while loop (don't use the
 *  elusive do-while loop, the semantics of calling Done after an iteration
 *  F things up).  You should only use which ever loop fits
 *  your current parallel need best, and not both loops.
 *
 *  \code
 *  //Assuming a continuation of the code above
 *  MPIJob<int> Job(TaskQue);
 *
 *  //Example using Job in a for loop
 *  for(int i=Job.Begin();!Job.Done();i=Job.Next()){
 *      //Do the stuff you need to for task i
 *  }
 *
 *
 *  //Example using Job in a while loop
 *
 *  int i=Job.Begin();
 *  while(!Job.Done()){
 *     //Do the stuff for task i
 *
 *     i=Job.Next();
 *  }
 *  \endcode
 *
 *
 *  After an MPI process completes its loop and calls Done() for the final
 *  time (i.e. the time Done() returns true), that process is free to continue
 *  on and do whatever it wants.  I have modified the WorldComm object
 *  so that it doesn't forget the communicators, so as long as you don't
 *  let your MPIJob object go out of scope you can still synchronize
 *  your data down the road.  This provides you the opportunity to start
 *  another task while you wait for the former to end.
 *
 *  Once you are ready to get your data back or do something like an all
 *  reduce (add up all the distributed numbers), you now call the approrpiate
 *  MPIJob function to get your data.  Look at the available calls to see
 *  what's available.
 *
 *
 *  A basically complete example for our interaction energy:
 *
 *  \code
 *  //I'm lazy so still using typedefs
 *  typedef MPITask<std::string> TaskType;
 *  typedef boost::shared_ptr<TaskType> SharedTask;
 *
 *  //Make our tasks, obviously you want a fxn for more complicated cases
 *  std::vector<SharedTask> Tasks;
 *  Tasks.push_back(SharedTask(new TaskType("Complex",0));
 *  Tasks.push_back(SharedTask(new TaskType("Chloride",1));
 *  Tasks.push_back(SharedTask(new TaskType("Water",2));
 *
 *  //Set-up MPI
 *  MPIJob<std::string> Job(Tasks);
 *
 *  //Calculate our energies
 *  std::vector<double> Energies;
 *  for(std::string job=Job.Begin();!Job.Done();job=Job.Next()){
 *      //We will pretend this function returns the appropriate molecule
 *      //given it's ASCII name
 *      molecule mol=Set_Molecule(job);
 *
 *      //We will pretend that this calculates the energy
 *      double energy=CalcEnergy(mol);
 *
 *      //Now we save our energy to our local array
 *      Energies.push_back(energy);
 *
 *  }
 *  //Call synch to make sure every process has all values
 *  std::vector<double> AllEnergies=Job.Synch(Energies,1);
 *
 *  //And finally actually compute our energy
 *  std::cout<<"The interaction energy is: "<<
 *              AllEnergies[0]-AllEnergies[1]-AllEnergies[2]<<std::endl;
 *
 *  \endcode
 *
 *  \param[in] T The type (int,string, etc.) of the user's label)
 *
 */
template <typename T>
class MPIJob:protected LibParallel::MPIJobGuts {
   private:
      typedef MPITask<T> TaskType;
      typedef LibParallel::MPIJobGuts BaseType;
   protected:
      ///The tasks the user gave us
      std::vector<TaskType> Tasks_;

      ///Default constructor, protected so only derived classes may call it
      MPIJob() {
      }

      /** \brief Constructor that will be called by the Advanced version of
       *         this class.
       *
       *         This function is essentially a backdoor to the base class,
       *         for the derived class.  In reality this should probably be
       *         done via multiple inheritance, but I'm trying to avoid that.
       *         For all intents and purposes, this function is irrelevant
       *         for the current class.
       */
      MPIJob(const std::vector<TaskType>& Tasks, int MaxProcs,
            bool ForceDynamic) :
            MPIJobGuts(Tasks, MaxProcs, ForceDynamic), Tasks_(Tasks) {
      }

   public:
      ///No memory to free-up
      virtual ~MPIJob<T>() {
      }

      /** \brief Creates an MPIJob object
       *
       * Given a set of tasks to perform, picks the best MPI algorithm
       * based on the priorities of the jobs.
       *
       * \param[in] Tasks An array of your desired tasks
       */
      MPIJob<T>(const std::vector<TaskType>& Tasks) :
            MPIJobGuts(Tasks), Tasks_(Tasks) {
         if (Tasks_.size()<1)Error("You can't parallelize 0 tasks!!!");
      }

      /** \brief Returns the next job you should run
       *
       *   In particular note that if you haven't started iterating a call
       *   to next is equivalent to a call to Begin; however, they are not
       *   analogous in any other circumstance.
       */
      T Next() {
         return Tasks_[BaseType::Next()].GetLabel();
      }

      /** \brief Returns true if this MPI process is done
       *
       * If you have no tasks to run, this will be true from the get
       * go.
       */
      bool Done() {
         return BaseType::Done();
      }

      /** \brief Returns the first job you should run
       *
       *   One thing to note about the begin function is that calling
       *   it twice, or after calling Next(), will re-add all previously
       *   completed jobs to your que.  Think of this analogous to looping
       *   over a std::vector using an iterator set to the beginning of the
       *   vector with std::vector::begin().  If you call it twice, you
       *   start looping all over again.
       */
      T Begin() {
         return Tasks_[BaseType::Begin()].GetLabel();
      }

      /** \brief Takes distributed values and broadcasts them to all
       *         processes.
       *
       *  This function's description is easier if we start with:
       *
       *  \param[in] N The number of elements in LocalValues that come from
       *               each task
       *
       *  If our current MPI process ran n tasks, then LocalValues is a
       *  n*N long array.  The first N came from the first task you were
       *  given, the second N from the second task you were given, etc.
       *
       *  \param[in] LocalValues The chunk of the distributed data this MPI
       *             Task collected.
       *
       *  Given this data, MPIJob will put together an array and return it.
       *  This array will be N*TotalNumberOfTasks long, and the order will
       *  be whatever order your MPITasks were in when you initialized this
       *  MPIJob object
       *
       *  \return A local copy of all elements, arranged according to input
       *          order
       *
       */
      template <typename T2>
      std::vector<T2> Synch(const std::vector<T2>& LocalValues, const int N) {
         return BaseType::Synch(LocalValues, N);
      }

      /** \brief Takes a vector of values and combines them element wise
       *         across MPI processes.
       *
       *  The explanation for how to use this, probably works best by
       *  example.
       *
       *  Say you are parallelizing the determination of the center of mass
       *  of the molecule, the formula for which is:
       *  \f[
       *    r_i=\frac{\sum_j m_jr_{ji}}{\sum_j m_j},
       *  \f]
       *  where \f$r_i\f$ is the i-th component of the center of mass,
       *  \f$m_j\f$ is the mass of the j-th atom, and \f$r_{ji}\f$ is the
       *  position of the j-th atom, in the i-th direction.
       *
       *  To parallelize it you distribute the atoms over MPI processes,
       *  then on each MPI process you calculate four things, the numerator
       *  for each component and the denominator, for each atom you are in
       *  charge of.  You then have two reduce calls:
       *  \code
       *  std::vector<double> Numerators=reduce(LocalNumerators,3,ADD);
       *  std::vector<double> Denominator=reduce(LocalDenominators,1,ADD);
       *  for(int i=0;i<3;i++)
       *     CoM[i]=Numerators[i]/Denominator;
       *  \endcode
       *
       *  where hopefully the variables are self-explanatory.
       *
       *  \param[in] LocalValues Your local contribution to the final value
       *  \param[in] N The number of elements in the final result
       *  \param[in] op How the results are being combined, can be any basic
       *                binary operation: +,-,*,/,%. See parallel.h for
       *                all supported operations.
       *
       */
      template <typename T2>
      std::vector<T2> Reduce(const std::vector<T2>& LocalValues, const int N,
            const LibParallel::MPIOperation& op) const {
         return BaseType::Reduce(LocalValues, N, op);
      }

};

/** \brief A more advanced MPIJob class, this entire class is considered
 *         expert.
 *
 * This class can be used as a wrapper to MPIJob, simply by providing
 * no arguments to the constructor.
 *
 * The first optional parameter is the maximum number of MPI processes
 * this task this Job may use.  A negative value says you want all
 * processes available to you, 0 will prevent the job from getting any
 * MPI processes (don't know why you'd want to do this, but the code
 * will allow it, maybe as a quick way of turning of parallel regions?),
 * and a positive, non-zero, value sets the maximum upper-bound your job
 * will use (MaxProcs=2 gives you at most 2 processes).  It is
 * important to note that you may get less processes
 * than the value of MaxProcs; however, you will never get more. What
 * makes this an advanced option is that you need to be more involved
 * in the MPI process, and this class loses some automagic.
 *
 * For the interested reader, let me elaborate on the last sentence.
 * Let's assume we are starting with Comm1, which has a total of 6
 * processes (numbered 0 through 5 on Comm1), and we only want to use
 * 3 processes.  We make Comm2, by placing 0,1, and 2
 * in one communicator, and 3,4, and 5 in another.  Our Comms look
 * like:
 *
 * \code{.unparsed}
 *  Comm1:     0-5
 *              |
 *             / \
       *            /   \
       *           /     \
       * Comm2:  0-2     3-5
 * \endcode{}
 *
 * Processes 3-5 will "fall through" your loop by setting their
 * Done() routine to true.
 *
 * \code
 * MPIJobAdvanced<T> Job1(TaskSet1,3);
 * if(Job1.EnoughProcs()){
 *    if(!Job1.Done()){
 *        for(int task=Job1.Begin();!Job1.Done();task=Job1.Next()){
 *            //Do stuff that only needs 3 processors
 *        }
 *    }
 *    else{
 *        //When this gets created it will be on Comm2 with processes 3-5
 *        //If we wanted this could be another MPIJobAdvanced and we
 *        //could further split 3-5
 *        MPIJob Job2(TaskSet2);
 *        for(int task=Job2.Begin();!Job2.Done();task=Job2.Next(){
 *            //Do stuff that uses rest of processors
 *        }
 *    }
 * }
 * else{//You better have a back-up plan for scenarios where you don't
 *      //have enough processors available
 * }
 * \endcode
 *
 * Note, that if we had 3 or fewer processors the else clause of the
 * if(!Job1.Done()) statement would never have been reached because
 * all processors would have gone into the first half of the if.  This
 * is the point of the outermost if statement, to catch cases where
 * the else wouldn't have triggered.
 *
 */
template <typename T>
class MPIJobAdvanced:public MPIJob<T> {
   private:
      typedef MPITask<T> TaskType;
   public:
      /** \brief Allows for a user to dynamically add tasks
       *
       *  Under dynamic task assignment models, such as Master/Slave, it
       *  is possible to have the scheduling MPI process add more tasks
       *  to the task queue before it is depleted, and while processes are
       *  still running.  Calling this function under any static scheduler
       *  will result in a PSIEXCEPTION.  To ensure that a dynamic scheduler
       *  is being used the user may explicitly one.  See the constructor
       *  for more details.
       */
      void AddTasks(const std::vector<TaskType>& Tasks) {
      }

      /** \brief A class for users who know some MPI and want to do more
       *        advanced things
       *
       *
       *  \param[in] Tasks The task we are performing
       *  \param[in] MaxProcs Maximum number of MPI processes to use
       *            Negative value means as many as are available
       *            (default is -1)
       *  \param[in] ForceDynamic If true a dynamic scheduling algorithm is
       *            necessarily used. False allows the MPIJobGuts to
       *            determine if a dynamic or a static scheduler is optimum.
       */
      MPIJobAdvanced<T>(const std::vector<TaskType>& Tasks, const int MaxProcs=
            -1, const bool ForceDynamic=false) :
            MPIJob<T>(Tasks, MaxProcs, ForceDynamic) {
                LibParallel::MPIJobGuts::Enough_=LibParallel::MPIJobGuts::EnoughProcs();
      }

      ///Returns true if there are enough processes for your MaxProc request
      bool EnoughProcs() const {
         return LibParallel::MPIJobGuts::EnoughProcs();
      }

      ///Closes off the area where the processes are allowed to seperate
      void Wait() const {
         LibParallel::MPIJobGuts::Wait();
      }
};

}

#endif /* MPIJOB_H_ */
