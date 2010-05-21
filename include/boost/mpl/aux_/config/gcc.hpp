
#ifndef BOOST_MPL_AUX_CONFIG_GCC_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_GCC_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: gcc.hpp 49239 2008-10-10 09:10:26Z agurtovoy $
// $Date: 2008-10-10 05:10:26 -0400 (Fri, 10 Oct 2008) $
// $Revision: 49239 $

#if defined(__GNUC__) && !defined(__EDG_VERSION__)
#   define BOOST_MPL_CFG_GCC ((__GNUC__ << 8) | __GNUC_MINOR__)
#else
#   define BOOST_MPL_CFG_GCC 0
#endif

#endif // BOOST_MPL_AUX_CONFIG_GCC_HPP_INCLUDED
