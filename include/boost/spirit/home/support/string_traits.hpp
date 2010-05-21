/*=============================================================================
    Copyright (c) 2001-2010 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
================================================_==============================*/
#if !defined(BOOST_SPIRIT_STRING_TRAITS_OCTOBER_2008_1252PM)
#define BOOST_SPIRIT_STRING_TRAITS_OCTOBER_2008_1252PM

#if defined(_MSC_VER)
#pragma once
#endif

#include <boost/spirit/home/support/container.hpp>
#include <string>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/utility/enable_if.hpp>

namespace boost { namespace spirit { namespace traits
{
    ///////////////////////////////////////////////////////////////////////////
    // Determine if T is a string
    ///////////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_string : mpl::false_ {};

    template <typename T>
    struct is_string<T const> : is_string<T> {};

    template <>
    struct is_string<char const*> : mpl::true_ {};

    template <>
    struct is_string<wchar_t const*> : mpl::true_ {};

    template <>
    struct is_string<char*> : mpl::true_ {};

    template <>
    struct is_string<wchar_t*> : mpl::true_ {};

    template <size_t N>
    struct is_string<char[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<wchar_t[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<char const[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<wchar_t const[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<char(&)[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<wchar_t(&)[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<char const(&)[N]> : mpl::true_ {};

    template <size_t N>
    struct is_string<wchar_t const(&)[N]> : mpl::true_ {};

    template <typename T, typename Traits, typename Allocator>
    struct is_string<std::basic_string<T, Traits, Allocator> > : mpl::true_ {};

    ///////////////////////////////////////////////////////////////////////////
    // Get the underlying char type of a string
    ///////////////////////////////////////////////////////////////////////////
    template <typename T>
    struct char_type_of;

    template <typename T>
    struct char_type_of<T const> : char_type_of<T> {};

    template <>
    struct char_type_of<char const*> : mpl::identity<char const> {};

    template <>
    struct char_type_of<wchar_t const*> : mpl::identity<wchar_t const> {};

    template <>
    struct char_type_of<char*> : mpl::identity<char> {};

    template <>
    struct char_type_of<wchar_t*> : mpl::identity<wchar_t> {};

    template <size_t N>
    struct char_type_of<char[N]> : mpl::identity<char> {};

    template <size_t N>
    struct char_type_of<wchar_t[N]> : mpl::identity<wchar_t> {};

    template <size_t N>
    struct char_type_of<char const[N]> : mpl::identity<char const> {};

    template <size_t N>
    struct char_type_of<wchar_t const[N]> : mpl::identity<wchar_t const> {};

    template <size_t N>
    struct char_type_of<char(&)[N]> : mpl::identity<char> {};

    template <size_t N>
    struct char_type_of<wchar_t(&)[N]> : mpl::identity<wchar_t> {};

    template <size_t N>
    struct char_type_of<char const(&)[N]> : mpl::identity<char const> {};

    template <size_t N>
    struct char_type_of<wchar_t const(&)[N]> : mpl::identity<wchar_t const> {};

    template <typename T, typename Traits, typename Allocator>
    struct char_type_of<std::basic_string<T, Traits, Allocator> >
      : mpl::identity<T> {};

    ///////////////////////////////////////////////////////////////////////////
    // Get the C string from a string
    ///////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline T* get_c_string(T* str) { return str; }

    template <typename T>
    inline T const* get_c_string(T const* str) { return str; }

    template <typename T, typename Traits, typename Allocator>
    inline T const* get_c_string(std::basic_string<T, Traits, Allocator>& str)
    { return str.c_str(); }

    template <typename T, typename Traits, typename Allocator>
    inline T const* get_c_string(std::basic_string<T, Traits, Allocator> const& str)
    { return str.c_str(); }

    ///////////////////////////////////////////////////////////////////////////
    // Get the begin/end iterators from a string
    ///////////////////////////////////////////////////////////////////////////

    // Implementation for C-style strings.
    template <typename T>
    inline T const* get_begin(T const* str) { return str; }

    template <typename T>
    inline T* get_begin(T* str) { return str; }

    template <typename T>
    inline T const* get_end(T const* str)
    {
        T const* last = str;
        while (*last)
            last++;
        return last;
    }

    template <typename T>
    inline T* get_end(T* str)
    {
        T* last = str;
        while (*last)
            last++;
        return last;
    }

    // Implementation for containers (includes basic_string).
    template <typename T, typename Str>
    inline typename Str::const_iterator get_begin(Str const& str)
    { return str.begin(); }

    template <typename T, typename Str>
    inline typename Str::iterator get_begin(Str& str)
    { return str.begin(); }

    template <typename T, typename Str>
    inline typename Str::const_iterator get_end(Str const& str)
    { return str.end(); }

    template <typename T, typename Str>
    inline typename Str::iterator get_end(Str& str)
    { return str.end(); }

    // Default implementation for other types: try a C-style string
    // conversion.
    // These overloads are explicitly disabled for containers,
    // as they would be ambiguous with the previous ones.
    template <typename T, typename Str>
    inline typename disable_if<is_container<Str>
      , T const*>::type get_begin(Str const& str)
    { return str; }

    template <typename T, typename Str>
    inline typename disable_if<is_container<Str>
      , T const*>::type get_end(Str const& str)
    { return get_end(get_begin<T>(str)); }
}}}

#endif
