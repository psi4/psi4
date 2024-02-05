/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_bin_psimrcc_idmrpt2_h
#define _psi_src_bin_psimrcc_idmrpt2_h
/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "blas.h"
#include "manybody.h"
#include "psimrcc_wfn.h"
#include "updater.h"
#include "psi4/liboptions/liboptions.h"

namespace psi {
namespace psimrcc {

class Updater;

/**
        @author Francesco Evangelista <frank@ccc.uga.edu>
*/
class IDMRPT2 : public CCManyBody {
   public:
    IDMRPT2(std::shared_ptr<PSIMRCCWfn> wfn, Options& options);
    ~IDMRPT2() override;
    double compute_energy() override;

   private:
    std::shared_ptr<Updater> updater_;

    void add_matrices();
    void read_mrpt2_integrals();
    void update_amps_mkpt2();
    void synchronize_amps();
    void build_amplitudes();
    void build_t1_ia_amplitudes();
    void build_t1_IA_amplitudes();
    void build_t2_ijab_amplitudes();
    void build_t2_iJaB_amplitudes();
    void build_t2_IJAB_amplitudes();

    void build_Heff_mrpt2_diagonal();
    void build_Heff_scs_mrpt2_diagonal();
    void build_Heff_mrpt2_offdiagonal();
    void build_Heff_uv();
    void build_Heff_UV();
    void build_Heff_uVxY();
    void build_Heff_uvxy();
    void build_Heff_UVXY();

    void build_F_intermediates();
    void build_F_ae_intermediates();
    void build_F_AE_intermediates();
    void build_F_mi_intermediates();
    void build_F_MI_intermediates();
};

}  // namespace psimrcc
}  // namespace psi

#endif  // _psi_src_bin_psimrcc_idmrpt2_h
