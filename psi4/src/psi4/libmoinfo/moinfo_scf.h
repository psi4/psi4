/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmoinfo_moinfo_scf_h_
#define _psi_src_lib_libmoinfo_moinfo_scf_h_

/*! \file
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#include <string>

#include "psi4/libmints/wavefunction.h"
#include "moinfo_base.h"

namespace psi {

class MOInfoSCF : public MOInfoBase {
   public:
    MOInfoSCF(Wavefunction& ref_wfn_, Options& options_, bool silent_ = false);
    ~MOInfoSCF();

    bool get_guess_occupation_flag() const { return (guess_occupation_flag); }

   private:
    void read_mo_spaces();
    void print_mo();
    bool guess_occupation_flag;
};

extern MOInfoSCF* moinfo_scf;

}  // namespace psi

#endif  // _psi_src_lib_libmoinfo_moinfo_scf_h_
