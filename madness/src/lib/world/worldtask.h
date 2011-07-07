/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: worldtask.h 2395 2011-07-02 02:36:51Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDTASK_H__INCLUDED
#define MADNESS_WORLD_WORLDTASK_H__INCLUDED

/// \file worldtask.h
/// \brief Defines TaskInterface and implements WorldTaskQueue and associated stuff.

// To do:
// a) Redo this ontop of serializable tasks which will remote much of the clutter
//    due to multiple length argument lists.
// b) Stealing which pretty much presume a) has been done

#include <iostream>
#include <world/nodefaults.h>
//#include <world/worldtypes.h>
//#include <world/typestuff.h>
#include <world/worlddep.h>
#include <world/worldthread.h>
#include <world/worldrange.h>
#include <world/worldfut.h>
#include <world/worldtime.h>

namespace madness {

    // Forward decls
    class World;
    class WorldTaskQueue;
    template <typename functionT> struct TaskFunction;
    template <typename memfunT> struct TaskMemfun;

    /// All world tasks must be derived from this public interface

    /// Multiple worlds with independent queues feed tasks into shared task
    /// pool that is mapped to the H/W.
    ///
    /// For simplicity and backward compatibility we maintain two run interfaces
    /// but new code should adopt the multithreaded interface
    ///
    /// \c run(World&) - the user implements this for a single-threaded task
    ///
    /// \c run(World&, \c const \c TaskThreadEnv&) - the user implements this for
    /// a multi-threaded task.
    ///
    class TaskInterface : public DependencyInterface , public PoolTaskInterface {
        friend class WorldTaskQueue;
    private:
        volatile World* world;
        CallbackInterface* completion;

        // Used for submission to underlying queue when all dependencies are satisfied
        struct Submit : public CallbackInterface {
            PoolTaskInterface* p;
            Submit(PoolTaskInterface* p) : p(p) {}
            void notify() {
                ThreadPool::add(p);
            }
        } submit;


        /// Set task info

        /// \param w The world object that contains the task
        /// \param c Call this callback on completion
        void set_info(World* w, CallbackInterface* c) {
            world = w;
            completion = c;
        }

        /// Adds call back to schedule task when outstanding dependencies are satisfied
        void register_submit_callback() { register_callback(&submit); }

    protected:
        virtual void run(const TaskThreadEnv& env);

    public:
        static bool debug;

        /// Create a new task with ndepend dependencies (default 0) and given attributes
        TaskInterface(int ndepend=0, const TaskAttributes attr = TaskAttributes())
                : DependencyInterface(ndepend)
                , PoolTaskInterface(attr)
                , world(0)
                , completion(0)
                , submit(this)
        {}

        /// Create a new task with zero dependencies and given attributes
        explicit TaskInterface(const TaskAttributes& attr)
                : DependencyInterface(0)
                , PoolTaskInterface(attr)
                , world(0)
                , completion(0)
                , submit(this)
        {}

//         void serialize(Buffer& ar) {
//             throw "there is no way this is correct";
//             ar & *static_cast<PoolTaskInterface*>(this) & world;
//         }

        /// Runs a single-threaded task ... derived classes must implement this.

        /// This interface may disappear so new code should use the multi-threaded interface.
        virtual void run(World&) {
            //print("in virtual run(world) method");
            MADNESS_EXCEPTION("World TaskInterface: user did not implement one of run(world) or run(world, taskthreadenv)", 0);
        }

        /// Runs a multi-threaded task
        virtual void run(World& world, const TaskThreadEnv& env) {
            //print("in virtual run(world,env) method", env.nthread(), env.id());
            if (env.nthread() != 1)
                MADNESS_EXCEPTION("World TaskInterface: user did not implement run(world, taskthreadenv) for multithreaded task", 0);
            run(world);
        }

        World* get_world() const { return const_cast<World*>(world); }

        virtual ~TaskInterface() { if (completion) completion->notify(); }

    };

#include <world/taskfn.h>

    /// Multi-threaded queue to manage and run tasks.
    class WorldTaskQueue : public CallbackInterface, private NO_DEFAULTS {
        friend class TaskInterface;
    private:
        World& world;              ///< The communication context
        const ProcessID me;        ///< This process
        AtomicInt nregistered;     ///< Counts pending tasks

        void notify() { nregistered--; }

        // Used in for_each kernel to check completion
        static bool completion_status(bool left, bool right) {
            return (left && right);
        }


        // Used in reduce kernel
        template <typename resultT, typename opT>
        static resultT sum(const resultT& left, const resultT& right, const opT& op) {
            //std::cout << " REDUCE SUM " << left << " " << right << std::endl;
            return op(left,right);
        }


    public:
        WorldTaskQueue(World& world);

        /// Returns the number of pending tasks
        size_t size() const { return nregistered; }

        /// Returns the number of threads in the pool
        size_t nthread() { return ThreadPool::size(); }

        /// Returns all of the threads in the thread pool
        pthread_t* thread_id() { return ThreadPool::thread_id(); }

        /// Returns the thread id of the communicator thread
        pthread_t comm_thread_id() { return RMI::get_comm_id(); }

        /// Add a new local task taking ownership of the pointer

        /// The task pointer (t) is assumed to have been created with
        /// \c new and when the task is eventually run the queue
        /// will call the task's destructor using \c delete.
        ///
        /// Once the task is complete it will execute
        /// task_complete_callback to decrement the number of pending
        /// tasks and be deleted.
        void add(TaskInterface* t)  {
            nregistered++;

            t->set_info(&world, this);       // Stuff info

            if (t->ndep() == 0) {
                ThreadPool::add(t); // If no dependencies directly submit
            } else {
                // With dependencies must use the callback to avoid race condition
                t->register_submit_callback();
                //t->dec();
            }
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T,
            typename a4T, typename a5T, typename a6T, typename a7T, typename a8T,
            typename a9T>
        typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::futureT
        add(TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>* t) {

            typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::futureT
                res(t->result());
            add(static_cast<TaskInterface*>(t));
            return res;
        }

        /// Reduce op(item) for all items in range using op(sum,op(item))

        /// The operation must provide the following interface of
        /// which the \c operator() methods are required by reduce()
        /// and the rest by the task interface.
        /// \code
        /// struct opT {
        ///     opT();
        ///     opT(const &opT);
        ///     resultT operator()(const rangeT::iterator& it) const;
        ///     resultT operator()(const resultT& left, const resultT& right);
        ///     template <typename Archive> void serialize(const Archive& ar);
        /// }
        /// \endcode
        /// Note that the serialize method does not actually have to
        /// work unless you want to have the task be stealable.
        /// Adjust the chunksize in the range to control granularity.
        template <typename resultT, typename rangeT, typename opT>
        Future<resultT> reduce(const rangeT& range, const opT& op) {
            if (range.size() <= range.get_chunksize()) {
                resultT sum = resultT();
                for (typename rangeT::iterator it=range.begin(); it != range.end(); ++it) sum = op(sum,op(it));
                return Future<resultT>(sum);
            } else {
                rangeT left = range;
                rangeT right(left,Split());

                Future<resultT>  leftsum = add(*this, &WorldTaskQueue::reduce<resultT,rangeT,opT>, left,  op);
                Future<resultT> rightsum = add(*this, &WorldTaskQueue::reduce<resultT,rangeT,opT>, right, op);
                return add(&WorldTaskQueue::sum<resultT,opT>, leftsum, rightsum, op);
            }
        }

        /// Apply op(item) for all items in range

        /// The operation must provide the following interface of
        /// which the \c operator() method is required by for_each()
        /// and the rest by the task interface.
        /// \code
        /// struct opT {
        ///     opT();
        ///     opT(const &opT);
        ///     bool operator()(const rangeT::iterator& it) const;
        ///     template <typename Archive> void serialize(const Archive& ar);
        /// }
        /// \endcode
        /// Note that the serialize method does not actually have to
        /// work unless you want to have the task be stealable.
        ///
        /// Adjust the chunksize in the range to control granularity.
        ///
        /// Your operation should return true/false for success failure
        /// and the logical and of all results is returned as the
        /// future result.
        ///
        /// You can ignore the result if you are interested
        /// in neither synchronization nor result status.
        template <typename rangeT, typename opT>
        Future<bool> for_each(const rangeT& range, const opT& op) {
            if (range.size() <= range.get_chunksize()) {
                bool status = true;
                for (typename rangeT::iterator it=range.begin();  it != range.end();  ++it) status &= op(it);
                return Future<bool>(status);
            } else {
                rangeT left = range;
                rangeT right(left,Split());
                Future<bool>  leftsum = add(*this, &WorldTaskQueue::for_each<rangeT,opT>, left,  op);
                Future<bool> rightsum = add(*this, &WorldTaskQueue::for_each<rangeT,opT>, right, op);
                return add(&WorldTaskQueue::completion_status, leftsum, rightsum);
            }
        }


        /// Invoke "resultT (*function)(void)" as a local task

        /// A future is returned to hold the eventual result of the task.
        /// Future<void> is an empty class that may be ignored.
        template <typename functionT>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(functionT function, const TaskAttributes attr=TaskAttributes()) {
            return add(me,function,attr);
        }


        /// Invoke "resultT (*function)(void)" as a task, local or remote

        /// A future is returned to hold the eventual result of the task.
        /// Future<void> is an empty class that may be ignored.
        template <typename functionT>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me)
                add(new TaskFunction<functionT>(result, function, attr));
            else
                TaskFunction<functionT>::sender(world, where, result, function, attr);
            return result;
        }


        /// Invoke "resultT (*function)(arg1T)" as a local task
        template <typename functionT, typename arg1T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(functionT function, const arg1T& arg1, const TaskAttributes attr=TaskAttributes()) {
            return add(me, function, arg1, attr);
        }


        /// Invoke "resultT (*function)(arg1T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        ///
        /// Arguments must be (de)serializable and must of course make
        /// sense at the remote destination.  Fundamental types,
        /// simple STL containers, and pointers to World,
        /// WorldContainer, and user-defined types derived from
        /// WorldObject<> are automatically handled.  Anything else is
        /// your problem.
        ///
        /// An argument that is a future may be used to carry
        /// dependencies for local tasks.  An unready future cannot be
        /// used as an argument for a remote tasks --- i.e., remote
        /// tasks must be ready to execute (you can work around this
        /// by making a local task to submit the remote task once
        /// everything is ready).
        template <typename functionT, typename arg1T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, attr);
            }
            return result;
        }


        /// Invoke "resultT (*function)(arg1T,arg2T)" as a local task
        template <typename functionT, typename arg1T, typename arg2T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(functionT function,
            const arg1T& arg1, const arg2T& arg2, const TaskAttributes attr=TaskAttributes()) {
            return add(me,function,arg1,arg2,attr);
        }


        /// Invoke "resultT (*function)(arg1T,arg2T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, attr);
            }
            return result;
        }


        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T)" as a local task
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const TaskAttributes attr=TaskAttributes()) {
            return add(me,function,arg1,arg2,arg3,attr);
        }


        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, attr);
            }
            return result;
        }

        class Stealer {
            WorldTaskQueue& q;
            std::vector<TaskInterface*>& v;
            const int nsteal;

        public:
            Stealer(WorldTaskQueue& q, std::vector<TaskInterface*>& v, int nsteal)
                : q(q)
                , v(v)
                , nsteal(nsteal)
            {}

            bool operator()(PoolTaskInterface** pt);
        };


        /// Invoke locally or remotely to send tasks to process P
        std::vector<TaskInterface*> steal(int nsteal);


        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, attr);
            }
            return result;
        }



        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, attr);
            }
            return result;
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, attr);
            }
            return result;
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg7));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
            }
            return result;
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6,
	    const arg7T& arg7, const arg8T& arg8, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
	      add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg7));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg8));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, attr);
            }
            return result;
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))>
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6,
	    	const arg7T& arg7, const arg8T& arg8,
			const arg9T& arg9, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
	      add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg7));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg8));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg9));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, attr);
            }
            return result;
        }


        /// Invoke "resultT (obj.*memfun)()" as a local task
        template <typename memfunT>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj, memfunT memfun, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T)" as a local task
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,attr));
            return result;
        }



        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,arg6,attr));
            return result;
        }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7,attr));
            return result;
        }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj,
            memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const arg8T& arg8,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8, attr));
            return result;
        }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8,arg9)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))>
        add(MEMFUN_OBJT(memfunT)& obj, memfunT memfun,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5,
			const arg6T& arg6, const arg7T& arg7, const arg8T& arg8, const arg9T& arg9,
            const TaskAttributes attr=TaskAttributes()) {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9, attr));
            return result;
        }

        struct ProbeAllDone {
            WorldTaskQueue* tq;
            double start;
            ProbeAllDone(WorldTaskQueue* tq) : tq(tq),start(cpu_time()) {}
            bool operator()() const;
        };

        /// Returns after all local tasks have completed

        /// While waiting the calling thread will run tasks.
        void fence();
    };


    /// Serialization container for sending tasks to remote nodes

    /// This is for internal use only. You should not use this class directly.
    /// \tparam refT The remote reference type for task result future
    /// \tparam functionT The task function type
    template <typename refT, typename functionT>
    struct TaskHandlerInfo {
        refT ref;               ///< Remote reference for a task result future
        functionT func;         ///< A task function
        TaskAttributes attr;    ///< Task attributes

        /// Construct task info object

        /// \param ref Remote reference to the result future
        /// \param func The task function
        /// \param attr The task attrubutes
        TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr)
                : ref(ref), func(func),attr(attr) {}
        TaskHandlerInfo() {}

        /// Serialization of object

        /// \tparam Archive The serialization archive type
        /// \param ar The serialization archive
        template <typename Archive>
        void serialize(const Archive& ar) {
            serialize_internal<functionT>(ar);
        }

    private:

        /// Identify function and member function pointers

        /// \tparam T The function to identify
        template <typename fnT>
        struct is_func_ptr {
            static const bool value =
                (std::is_function<typename std::remove_pointer<fnT>::type >::value
                || std::is_member_function_pointer<fnT>::value);
        };

        /// Serialization for function pointers and member function pointers

        /// \tparam fnT The function type
        /// \tparam Archive The serialization archive type
        /// \param ar The serialization archive
        template <typename fnT, typename Archive>
        typename enable_if<is_func_ptr<fnT> >::type
        serialize_internal(const Archive& ar) {
            ar & ref & archive::wrap_opaque(func) & attr;
        }

        /// Serialization for non- function pointers and member function pointers.

        /// \tparam fnT The function type
        /// \tparam Archive The serialization archive type
        /// \param ar The serialization archive
        template <typename fnT, typename Archive>
        typename disable_if<is_func_ptr<fnT> >::type
        serialize_internal(const Archive& ar) {
            ar & ref & func & attr;
        }
    };

    // Internal: Common functionality for TaskFunction and TaskMemfun classes
    class TaskFunctionBase : public TaskInterface {
    public:

        TaskFunctionBase(const madness::TaskAttributes& attributes)
                : TaskInterface(attributes) {}

        // Register non-ready future as a dependency
        template <typename T>
        inline void check_dependency(Future<T>& fut) {
            if (!fut.probe()) {
                inc();
                fut.register_callback(this);
            }
        }

        virtual ~TaskFunctionBase() {}
    };



    // Task wrapping "resultT (*function)()"
    template <typename resultT>
    struct TaskFunction<resultT(*)()> : public TaskFunctionBase {
        typedef resultT(*functionT)();
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg & info;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,info.attr));
        }

        static void sender(World& world, ProcessID dest, Future<REMFUTURE(resultT)>& result, functionT func, const TaskAttributes& attr) {
            world.am.send(dest, handler, new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr)));
        }

        futureT result;
        const functionT func;
        TaskFunction(const futureT& result, functionT func, const TaskAttributes& attr)
                : TaskFunctionBase(attr)
                , result(result)
                , func(func) {};

        void run(World& /*world*/) {
            result.set(func());
        }

        virtual ~TaskFunction() {}
    };


    // Task wrapping "resultT (*function)(arg1)"
    template <typename resultT, typename arg1_type>
    struct TaskFunction<resultT(*)(arg1_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg & info & arg1;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,info.attr));
        }

        template <typename a1T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& a1, const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(a1)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;

        template <typename a1T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& /*world*/) {
            result.set(func(arg1));
        }
    };


    // Task wrapping "resultT (*function)(arg1,arg2)"
    template <typename resultT, typename arg1_type, typename arg2_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type);

        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg & info & arg1 & arg2;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,info.attr));
        }

        template <typename a1T, typename a2T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2));
        }
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg & info & arg1 & arg2 & arg3;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2,arg3));
        }
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg & info & arg1 & arg2 & arg3 & arg4;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2),
                                     static_cast<arg3T>(arg3), static_cast<arg4T>(arg4)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2,arg3,arg4));
        }
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg & info & arg1 & arg2 & arg3 & arg4 & arg5;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3),
                                     static_cast<arg4T>(arg4), static_cast<arg5T>(arg5)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                     const a5T& a5, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2,arg3,arg4,arg5));
        }
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6,
                           const TaskAttributes& attr) {
            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3),
                                     static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
                     const a6T& a6, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2,arg3,arg4,arg5,arg6));
        }
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
    typename arg5_type, typename arg6_type, typename arg7_type>
    struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg7T arg7;
            arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,arg7,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6,
                           const a7T& arg7, const TaskAttributes& attr) {

            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3),
                                     static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6),
                                     static_cast<arg7T>(arg7)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                     const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
        }

        void run(World& /*world*/) {
            result.set(func(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
        }
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
      typename arg5_type, typename arg6_type, typename arg7_type, typename arg8_type>
      struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
		typedef REMFUTURE(REMCONST(REMREF(arg8_type))) arg8T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg7T arg7;
	    arg8T arg8;
            arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7 & arg8;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6,
                           const a7T& arg7, const a8T& arg8, const TaskAttributes& attr) {

            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3),
                                     static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6),
                                     static_cast<arg7T>(arg7), static_cast<arg8T>(arg8)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;
	Future<arg8T> arg8;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                     const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
	    	check_dependency(arg8);
        }

        void run(World& /*world*/) {
		  result.set(func(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8));
        }
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
      typename arg5_type, typename arg6_type, typename arg7_type, typename arg8_type, typename arg9_type>
      struct TaskFunction<resultT(*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type,arg9_type)> : public TaskFunctionBase {
        typedef resultT(*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type,arg9_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
		typedef REMFUTURE(REMCONST(REMREF(arg8_type))) arg8T;
		typedef REMFUTURE(REMCONST(REMREF(arg9_type))) arg9T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        using PoolTaskInterface::run;

        static void handler(const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg7T arg7;
	    	arg8T arg8;
	    	arg9T arg9;
            arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7 & arg8 & arg9;
            arg.get_world()->taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,info.attr));
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6,
                           const a7T& arg7, const a8T& arg8, const a9T& arg9, const TaskAttributes& attr) {

            world.am.send(dest, TaskFunction<functionT>::handler,
                          new_am_arg(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr),
                                     static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3),
                                     static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6),
                                     static_cast<arg7T>(arg7), static_cast<arg8T>(arg8), static_cast<arg9T>(arg9)));
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;
		Future<arg8T> arg8;
		Future<arg9T> arg9;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                     const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8), arg9(a9) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
		    check_dependency(arg8);
		    check_dependency(arg9);
        }

        void run(World& /*world*/) {
		  result.set(func(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9));
        }
    };


    // Task wrapping "resultT (obj.*function)()"
    template <typename resultT, typename objT>
    struct TaskMemfun<resultT(objT::*)()> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)();
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;

        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun) {}

        void run(World& /*world*/) {
            result.set((obj.*memfun)());
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1)"
    template <typename resultT, typename objT, typename arg1_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;

        template <typename a1T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const a2T& a2, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
    typename arg5_type, typename arg6_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
    typename arg5_type, typename arg6_type, typename arg7_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
      typename arg5_type, typename arg6_type, typename arg7_type, typename arg8_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
		typedef REMFUTURE(REMCONST(REMREF(arg8_type))) arg8T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;
		Future<arg8T> arg8;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
		    check_dependency(arg8);
        }

        void run(World& /*world*/) {
		  result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type,
      typename arg5_type, typename arg6_type, typename arg7_type, typename arg8_type, typename arg9_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type,arg9_type)> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type,arg9_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
		typedef REMFUTURE(REMCONST(REMREF(arg8_type))) arg8T;
		typedef REMFUTURE(REMCONST(REMREF(arg9_type))) arg9T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;
		Future<arg8T> arg8;
		Future<arg9T> arg9;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8), arg9(a9) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
		    check_dependency(arg8);
		    check_dependency(arg9);
        }

        void run(World& /*world*/) {
		  result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9));
        }
    };


    //
    // Same as above but now const
    //

    // Task wrapping "resultT (obj.*function)() const"
    template <typename resultT, typename objT>
    struct TaskMemfun<resultT(objT::*)()const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)()const;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;

        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun) {};

        void run(World& /*world*/) {
            result.set((obj.*memfun)());
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1)const"
    template <typename resultT, typename objT, typename arg1_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;

        template <typename a1T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1));
        }
    };


    // Task wrapping "resultT (obj.*function)(arg1,arg2)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)const> : public TaskFunctionBase {
        typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const TaskAttributes& attr)
                : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& /*world*/) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type, typename arg7_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type)const> : public TaskFunctionBase {
	    typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
		typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
		Future<arg7T> arg7;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
		    check_dependency(arg7);
        }

        void run(World& /*world*/) {
		  result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
        }
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type, typename arg7_type, typename arg8_type>
    struct TaskMemfun<resultT(objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type)const> : public TaskFunctionBase {
	    typedef resultT(objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type,arg8_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
		typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
		typedef REMFUTURE(REMCONST(REMREF(arg8_type))) arg8T;
        typedef Future<REMFUTURE(resultT)> futureT;

        using PoolTaskInterface::run;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
		Future<arg7T> arg7;
		Future<arg8T> arg8;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        TaskMemfun(const futureT& result, objT& obj, memfunT memfun,
                   const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const TaskAttributes& attr)
	  : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
		    check_dependency(arg7);
		    check_dependency(arg8);
        }

        void run(World& /*world*/) {
		  result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8));
        }
    };

}


#endif // MADNESS_WORLD_WORLDTASK_H__INCLUDED
