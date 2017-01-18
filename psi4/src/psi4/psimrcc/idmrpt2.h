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

#ifndef _psi_src_bin_psimrcc_idmrpt2_h
#define _psi_src_bin_psimrcc_idmrpt2_h
/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "manybody.h"
#include "psi4/liboptions/liboptions.h"

namespace psi{ namespace psimrcc{

class Updater;

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class IDMRPT2 : public CCManyBody
{
public:
  IDMRPT2(SharedWavefunction ref_wfn, Options &options);
  virtual ~IDMRPT2();
  void compute_mrpt2_energy(Updater* updater);
private:
  void add_matrices();
  void read_mrpt2_integrals();
  void update_amps_mkpt2(Updater* updater);
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

  void build_Heff_ijkabc();
  void build_Heff_ijKabC();
  void build_Heff_iJKaBC();
  void build_Heff_IJKABC();

  void build_F_intermediates();
  void build_F_ae_intermediates();
  void build_F_AE_intermediates();
  void build_F_mi_intermediates();
  void build_F_MI_intermediates();
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_idmrpt2_h
