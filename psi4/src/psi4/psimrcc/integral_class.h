/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef _psi_src_bin_psimrcc_integral_class_h_
#define _psi_src_bin_psimrcc_integral_class_h_

/*! \file
    \ingroup (PSIMRCC)
    \brief   This class is used to specify a class of integrals
*/

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class IntegralClass{
public:
  ///////////////////////////////////////////////////////////////////////////////
  // Class Constructor and Destructor
  ///////////////////////////////////////////////////////////////////////////////
  IntegralClass();
  ~IntegralClass();

  ///////////////////////////////////////////////////////////////////////////////
  // Class Interface
  ///////////////////////////////////////////////////////////////////////////////
  // Print routines
private:
  ///////////////////////////////////////////////////////////////////////////////
  // Class private functions
  ///////////////////////////////////////////////////////////////////////////////
  void        startup();
  void        cleanup();
  ///////////////////////////////////////////////////////////////////////////////
  // Class data
  ///////////////////////////////////////////////////////////////////////////////
  // Type                           // Name
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_integral_class_h_