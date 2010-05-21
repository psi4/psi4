
#ifndef BOOST_MPL_SIZE_T_HPP_INCLUDED
#define BOOST_MPL_SIZE_T_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2000-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: size_t.hpp 49239 2008-10-10 09:10:26Z agurtovoy $
// $Date: 2008-10-10 05:10:26 -0400 (Fri, 10 Oct 2008) $
// $Revision: 49239 $

#include <boost/mpl/size_t_fwd.hpp>

#define AUX_WRAPPER_VALUE_TYPE std::size_t
#define AUX_WRAPPER_NAME size_t
#define AUX_WRAPPER_PARAMS(N) std::size_t N

#include <boost/mpl/aux_/integral_wrapper.hpp>

#endif // BOOST_MPL_SIZE_T_HPP_INCLUDED
