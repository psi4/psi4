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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccsd_F_intr()
{

    // defs
    SharedTensor2d K, T, U, Tau;

    // OO block
    // F_mi =  \sum_{Q} t_Q b_mi^Q
    FijA->gemv(true, bQijA, T1c, 1.0, 0.0);

    // F_mi +=  \sum_{Q,e} Tau"_ie^Q b_me^Q
    Tau = SharedTensor2d(new Tensor2d("Tau2pp (Q|IA)", nQ, naoccA, navirA));
    Tau->read(psio_, PSIF_DFOCC_AMPS);
    FijA->contract332(false, true, navirA, bQiaA, Tau, 1.0, 1.0);
    Tau.reset();

    // VV block
    // F_ae =  \sum_{Q} t_Q b_ae^Q
    FabA->gemv(true, bQabA, T1c, 1.0, 0.0);

    // F_ae -=  \sum_{Q,m} Tau'_ma^Q b_me^Q
    Tau = SharedTensor2d(new Tensor2d("Tau2p (Q|IA)", nQ, naoccA, navirA));
    Tau->read(psio_, PSIF_DFOCC_AMPS);
    FabA->contract(true, false, navirA, navirA, nQ * naoccA, Tau, bQiaA, -1.0, 1.0);
    Tau.reset();

    // OV block
    // F_me +=  \sum_{Q} t_Q b_me^Q
    FiaA->gemv(true, bQiaA, T1c, 1.0, 0.0);

    // F_me -=  \sum_{Q,n} t_nm^Q b_ne^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    FiaA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, -1.0, 1.0);
    T.reset();

    // Ft_mi = F_mi + 1/2 \sum_{e} t_i^e F_me
    FtijA->gemm(false, true, FiaA, t1A, 0.5, 0.0);
    FtijA->add(FijA);

    // Ft_ae = F_ae - 1/2 \sum_{m} t_m^a F_me
    FtabA->gemm(true, false, t1A, FiaA, -0.5, 0.0);
    FtabA->add(FabA);

    //outfile->Printf("\tF int done.\n");

}// end ccsd_F_intr
}} // End Namespaces
