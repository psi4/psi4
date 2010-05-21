
#ifndef BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED
#define BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2001-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: arg_typedef.hpp 49239 2008-10-10 09:10:26Z agurtovoy $
// $Date: 2008-10-10 05:10:26 -0400 (Fri, 10 Oct 2008) $
// $Revision: 49239 $

#include <boost/mpl/aux_/config/lambda.hpp>
#include <boost/mpl/aux_/config/workaround.hpp>

#if defined(BOOST_MPL_CFG_NO_FULL_LAMBDA_SUPPORT) \
    || BOOST_WORKAROUND(__DMC__, BOOST_TESTED_AT(0x840))
    
#   define BOOST_MPL_AUX_ARG_TYPEDEF(T, name) typedef T name;

#else

#   define BOOST_MPL_AUX_ARG_TYPEDEF(T, name) /**/

#endif

#endif // BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED
