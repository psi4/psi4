#ifndef BOOST_THREAD_WIN32_ONCE_HPP
#define BOOST_THREAD_WIN32_ONCE_HPP

//  once.hpp
//
//  (C) Copyright 2005-7 Anthony Williams 
//  (C) Copyright 2005 John Maddock
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstring>
#include <cstddef>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/detail/interlocked.hpp>
#include <boost/thread/win32/thread_primitives.hpp>
#include <boost/thread/win32/interlocked_read.hpp>

#include <boost/config/abi_prefix.hpp>

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{
    using ::memcpy;
    using ::ptrdiff_t;
}
#endif

namespace boost
{
    struct once_flag
    {
        long status;
        long count;
        long throw_count;
        void* event_handle;

        ~once_flag()
        {
            if(count)
            {
                BOOST_ASSERT(count==throw_count);
            }
            
            void* const old_event=BOOST_INTERLOCKED_EXCHANGE_POINTER(&event_handle,0);
            if(old_event)
            {
                ::boost::detail::win32::CloseHandle(old_event);
            }
        }
    };

#define BOOST_ONCE_INIT {0,0,0,0}

    namespace detail
    {
        inline void* allocate_event_handle(void*& handle)
        {
            void* const new_handle=::boost::detail::win32::create_anonymous_event(
                ::boost::detail::win32::manual_reset_event,
                ::boost::detail::win32::event_initially_reset);
            
            void* event_handle=BOOST_INTERLOCKED_COMPARE_EXCHANGE_POINTER(&handle,
                                                                          new_handle,0);
            if(event_handle)
            {
                ::boost::detail::win32::CloseHandle(new_handle);
                return event_handle;
            }
            return new_handle;
        }
    }
    

    template<typename Function>
    void call_once(once_flag& flag,Function f)
    {
        // Try for a quick win: if the procedure has already been called
        // just skip through:
        long const function_complete_flag_value=0xc15730e2;
        long const running_value=0x7f0725e3;
        long status;
        bool counted=false;
        void* event_handle=0;
        long throw_count=0;

        while((status=::boost::detail::interlocked_read_acquire(&flag.status))
              !=function_complete_flag_value)
        {
            status=BOOST_INTERLOCKED_COMPARE_EXCHANGE(&flag.status,running_value,0);
            if(!status)
            {
                try
                {
                    if(!event_handle)
                    {
                        event_handle=::boost::detail::interlocked_read_acquire(&flag.event_handle);
                    }
                    if(event_handle)
                    {
                        ::boost::detail::win32::ResetEvent(event_handle);
                    }
                    f();
                    if(!counted)
                    {
                        BOOST_INTERLOCKED_INCREMENT(&flag.count);
                        counted=true;
                    }
                    BOOST_INTERLOCKED_EXCHANGE(&flag.status,function_complete_flag_value);
                    if(!event_handle && 
                       (::boost::detail::interlocked_read_acquire(&flag.count)>1))
                    {
                        event_handle=::boost::detail::allocate_event_handle(flag.event_handle);
                    }
                    if(event_handle)
                    {
                        ::boost::detail::win32::SetEvent(event_handle);
                    }
                    throw_count=::boost::detail::interlocked_read_acquire(&flag.throw_count);
                    break;
                }
                catch(...)
                {
                    if(counted)
                    {
                        BOOST_INTERLOCKED_INCREMENT(&flag.throw_count);
                    }
                    BOOST_INTERLOCKED_EXCHANGE(&flag.status,0);
                    if(!event_handle)
                    {
                        event_handle=::boost::detail::interlocked_read_acquire(&flag.event_handle);
                    }
                    if(event_handle)
                    {
                        ::boost::detail::win32::SetEvent(event_handle);
                    }
                    throw;
                }
            }

            if(!counted)
            {
                BOOST_INTERLOCKED_INCREMENT(&flag.count);
                counted=true;
                status=::boost::detail::interlocked_read_acquire(&flag.status);
                if(status==function_complete_flag_value)
                {
                    break;
                }
                event_handle=::boost::detail::interlocked_read_acquire(&flag.event_handle);
                if(!event_handle)
                {
                    event_handle=::boost::detail::allocate_event_handle(flag.event_handle);
                    continue;
                }
            }
            BOOST_VERIFY(!::boost::detail::win32::WaitForSingleObject(
                             event_handle,::boost::detail::win32::infinite));
        }
        if(counted || throw_count)
        {
            if(!BOOST_INTERLOCKED_EXCHANGE_ADD(&flag.count,(counted?-1:0)-throw_count))
            {
                if(!event_handle)
                {
                    event_handle=::boost::detail::interlocked_read_acquire(&flag.event_handle);
                }
                if(event_handle)
                {
                    BOOST_INTERLOCKED_EXCHANGE_POINTER(&flag.event_handle,0);
                    ::boost::detail::win32::CloseHandle(event_handle);
                }
            }
        }
    }
}

#include <boost/config/abi_suffix.hpp>

#endif
