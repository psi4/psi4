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

#ifndef thread_impl_h
#define thread_impl_h

#include "thread.h"
#include "runtime.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define DECLARE_STATIC_THREAD_WORKSPACE(x) \
static ThreadWorkspaceAllocatorTemplate<x> thr_workspace_##x; \
template<> x** ThreadWorkspaceAllocatorTemplate<x>::workspaces_ = 0


namespace yeti {

template <class T>
class ThreadWorkspaceAllocatorTemplate :
    public ThreadWorkspaceAllocator
{

    private:
        static T** workspaces_;

    public:
        ThreadWorkspaceAllocatorTemplate();

        static T* get_workspace(uli nthread);

        void allocate();

        void deallocate();

};

template <class T>
ThreadWorkspaceAllocatorTemplate<T>::ThreadWorkspaceAllocatorTemplate()
{
    ThreadEnvironment::add_allocator(this);
}

template <class T>
void
ThreadWorkspaceAllocatorTemplate<T>::allocate()
{
    uli nthread = YetiRuntime::nthread();
    workspaces_ = new T*[nthread];
    for (uli i=0; i < nthread; ++i)
    {
        workspaces_[i] = new T;
    }
}

template <class T>
void
ThreadWorkspaceAllocatorTemplate<T>::deallocate()
{
    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
        delete workspaces_[i];
    delete[] workspaces_;
}

template <class T>
T*
ThreadWorkspaceAllocatorTemplate<T>::get_workspace(uli threadnum)
{
    return workspaces_[threadnum];
}

template <class T>
T*
get_workspace(uli threadnum)
{
    return ThreadWorkspaceAllocatorTemplate<T>::get_workspace(threadnum);
}

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

