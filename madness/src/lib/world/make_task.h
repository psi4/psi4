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


 $Id $
 */

#ifndef MADNESS_WORLD_MAKE_TASK_H__INCLUDED
#define MADNESS_WORLD_MAKE_TASK_H__INCLUDED

#include <world/world.h>
#include <world/worldtask.h>

namespace madness {


    // Task function factories

    template <typename fnT>
    TaskFn<fnT>*
    make_task(fnT fn, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT> taskT;
        return new taskT(typename taskT::futureT(), fn, attr);
    }

    template <typename fnT, typename a1T>
    TaskFn<fnT, a1T>*
    make_task(fnT fn, const a1T& a1, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, attr);
    }

    template <typename fnT, typename a1T, typename a2T>
    TaskFn<fnT, a1T, a2T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T>
    TaskFn<fnT, a1T, a2T, a3T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
    TaskFn<fnT, a1T, a2T, a3T, a4T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
    TaskFn<fnT, a1T, a2T, a3T, a4T, a5T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, a5, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T>
    TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
            const a6T& a6, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, a5, a6, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T>
    TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
            const a6T& a6, const a7T& a7, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, a5, a6, a7, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T, typename a8T>
    TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
            const a6T& a6, const a7T& a7, const a8T& a8,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, a5, a6, a7, a8, attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T, typename a8T, typename a9T>
    TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>*
    make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
            const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;
        return new taskT(typename taskT::futureT(), fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr);
    }

    namespace {

        template <typename taskT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        void handler(const AmArg& arg) {
            MADNESS_ASSERT(taskT::arity <= 9u);

            // Get task info and arguments form active message

            TaskHandlerInfo<typename taskT::futureT::remote_refT,
                    typename taskT::functionT> info;
            a1T a1;
            a2T a2;
            a3T a3;
            a4T a4;
            a5T a5;
            a6T a6;
            a7T a7;
            a8T a8;
            a9T a9;

            arg & info & a1 & a2 & a3 & a4 & a5 & a6 & a7 & a8 & a9;

            // Create result future
            typename taskT::futureT result(info.ref);

            // Construct task
            taskT* task = NULL;
            switch(taskT::arity) {
            case 0u:
                task = new taskT(result, info.func, info.attr);
                break;
            case 1u:
                task = new taskT(result, info.func, a1, info.attr);
                break;
            case 2u:
                task = new taskT(result, info.func, a1, a2, info.attr);
                break;
            case 3u:
                task = new taskT(result, info.func, a1, a2, a3, info.attr);
                break;
            case 4u:
                task = new taskT(result, info.func, a1, a2, a3, a4, info.attr);
                break;
            case 5u:
                task = new taskT(result, info.func, a1, a2, a3, a4, a5, info.attr);
                break;
            case 6u:
                task = new taskT(result, info.func, a1, a2, a3, a4, a5, a6, info.attr);
                break;
            case 7u:
                task = new taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, info.attr);
                break;
            case 8u:
                task = new taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, a8, info.attr);
                break;
            case 9u:
                task = new taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, a8, a9, info.attr);
                break;
            }

            // Add task to queue
            arg.get_world()->taskq.add(task);
        }

        template <typename taskT, typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        inline typename taskT::futureT
        sender(World& world, ProcessID where, fnT fn, const a1T& a1,
                const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
                const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
                const TaskAttributes& attr)
        {
            typename taskT::futureT result;
            typedef TaskHandlerInfo<typename taskT::futureT::remote_refT, typename taskT::functionT> infoT;
            world.am.send(where, & handler<taskT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>,
                    new_am_arg(infoT(result.remote_ref(world), fn, attr),
                    a1, a2, a3, a4, a5, a6, a7, a8, a9));

            return result;
        }

        // Convert the task argument to the correct type for serialization.

        template <typename T>
        inline const T& am_arg(const Future<T>& f) {
            MADNESS_ASSERT(f.probe());
            return f.get();
        }

        template <typename T> inline const T& am_arg(const T& t) { return t; }

        typedef Future<void> voidT;

    } // namespace

    template <typename fnT>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, attr));

        return sender<taskT>(world, dest, fn, voidT(), voidT(), voidT(), voidT(),
                voidT(), voidT(), voidT(), voidT(), voidT(), attr);
    }

    template <typename fnT, typename a1T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), voidT(), voidT(),
                voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
    }

    template <typename fnT, typename a1T, typename a2T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn,  const a1T& a1, const a2T& a2,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), voidT(),
                voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn,  const a1T& a1, const a2T& a2,
            const a3T& a3, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, a5, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), am_arg(a5), voidT(), voidT(), voidT(), voidT(),
                attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, a5, a6, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), am_arg(a5), am_arg(a6), voidT(), voidT(), voidT(),
                attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, a5, a6, a7, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), voidT(), voidT(),
                attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, a5, a6, a7, a8, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                voidT(), attr);
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    typename task_result_type<fnT>::futureT
    send_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const a9T& a9, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;

        if(dest == world.rank())
            return world.taskq.add(make_task(fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr));

        return sender<taskT>(world, dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                am_arg(a9), attr);
    }


}  // namespace madness
#endif // MADNESS_WORLD_MAKE_TASK_H__INCLUDED
