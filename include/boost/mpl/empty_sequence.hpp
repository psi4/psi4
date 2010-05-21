
#ifndef BOOST_MPL_EMPTY_SEQUENCE_HPP_INCLUDED
#define BOOST_MPL_EMPTY_SEQUENCE_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2004
// Copyright Alexander Nasonov 2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: empty_sequence.hpp 49239 2008-10-10 09:10:26Z agurtovoy $
// $Date: 2008-10-10 05:10:26 -0400 (Fri, 10 Oct 2008) $
// $Revision: 49239 $

#include <boost/mpl/size_fwd.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/iterator_tags.hpp>

namespace boost { namespace mpl {

struct empty_sequence
{
    struct tag; 
    struct begin { typedef random_access_iterator_tag category; };    
    typedef begin end;
};

template<>
struct size_impl<empty_sequence::tag>
{
    template< typename Sequence > struct apply
        : int_<0>
    {
    };
};

}}

#endif // #ifndef BOOST_MPL_EMPTY_SEQUENCE_HPP_INCLUDED
