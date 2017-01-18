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

#ifndef _psi_src_bin_psimrcc_ccmatrixtmp_h
#define _psi_src_bin_psimrcc_ccmatrixtmp_h
/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/


namespace psi{ namespace psimrcc{

class CCMatrix;

enum DiskOpt {none,dump,release};

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCMatTmp{
public:
  CCMatTmp(CCMatrix* Matrix, DiskOpt disk_option = none);
  CCMatrix* operator->()   const {return(Matrix_);}
  CCMatrix* get_CCMatrix() const {return(Matrix_);}
  ~CCMatTmp();
private:
  DiskOpt     disk_option_;
  CCMatrix*   Matrix_;
};

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCMatIrTmp{
public:
  CCMatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option = none);
  CCMatrix* operator->()   const {return(Matrix_);}
  CCMatrix* get_CCMatrix() const {return(Matrix_);}
  ~CCMatIrTmp();
private:
  DiskOpt     disk_option_;
  int         irrep_;
  CCMatrix*   Matrix_;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccmatrixtmp_h