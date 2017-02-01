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

/** Standard library includes */
#include <fstream>
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::pair_index()
{

    // SO
    so_idx = SharedTensor2i(new Tensor2i("so_idx", nso_, nso_));
    for (int i = 0; i < nso_; i++) {
         for (int j = 0; j < nso_; j++) {
              so_idx->set(i, j, j + (i * nso_));
         }
    }

    // Active OO
    ij_idxAA = SharedTensor2i(new Tensor2i("IJ_idx", naoccA, naoccA));
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < naoccA; j++) {
              ij_idxAA->set(i, j, j + (i * naoccA));
         }
    }

    // All OO
    oo_idxAA = SharedTensor2i(new Tensor2i("OO_idx", noccA, noccA));
    for (int i = 0; i < noccA; i++) {
         for (int j = 0; j < noccA; j++) {
              oo_idxAA->set(i, j, j + (i * noccA));
         }
    }

    // Active OV
    ia_idxAA = SharedTensor2i(new Tensor2i("IA_idx", naoccA, navirA));
    for (int i = 0; i < naoccA; i++) {
         for (int a = 0; a < navirA; a++) {
              ia_idxAA->set(i, a, a + (i * navirA));
         }
    }

    // All OV
    ov_idxAA = SharedTensor2i(new Tensor2i("OV_idx", noccA, nvirA));
    for (int i = 0; i < noccA; i++) {
         for (int a = 0; a < nvirA; a++) {
              ov_idxAA->set(i, a, a + (i * nvirA));
         }
    }

    // Active VO
    ai_idxAA = SharedTensor2i(new Tensor2i("AI_idx", navirA, naoccA));
    for (int a = 0; a < navirA; a++) {
         for (int i = 0; i < naoccA; i++) {
              ai_idxAA->set(a, i, i + (a * naoccA));
         }
    }

    // All VO
    vo_idxAA = SharedTensor2i(new Tensor2i("VO_idx", nvirA, noccA));
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              vo_idxAA->set(a, i, i + (a * noccA));
         }
    }

    // Active VV
    ab_idxAA = SharedTensor2i(new Tensor2i("AB_idx", navirA, navirA));
    for (int a = 0; a < navirA; a++) {
         for (int b = 0; b < navirA; b++) {
              ab_idxAA->set(a, b, b + (a * navirA));
         }
    }

    // All VV
    vv_idxAA = SharedTensor2i(new Tensor2i("VV_idx", nvirA, nvirA));
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              vv_idxAA->set(a, b, b + (a * nvirA));
         }
    }


 // UNRESTRICTED
 if (reference_ == "UNRESTRICTED") {
    // Active Oo
    ij_idxAB = SharedTensor2i(new Tensor2i("Ij_idx", naoccA, naoccB));
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < naoccB; j++) {
              ij_idxAB->set(i, j, j + (i * naoccB));
         }
    }

    // Active oO
    ij_idxBA = SharedTensor2i(new Tensor2i("iJ_idx", naoccB, naoccA));
    for (int i = 0; i < naoccB; i++) {
         for (int j = 0; j < naoccA; j++) {
              ij_idxBA->set(i, j, j + (i * naoccA));
         }
    }

    // Active oo
    ij_idxBB = SharedTensor2i(new Tensor2i("ij_idx", naoccB, naoccB));
    for (int i = 0; i < naoccB; i++) {
         for (int j = 0; j < naoccB; j++) {
              ij_idxBB->set(i, j, j + (i * naoccB));
         }
    }

    // All Oo
    oo_idxAB = SharedTensor2i(new Tensor2i("Oo_idx", noccA, noccB));
    for (int i = 0; i < noccA; i++) {
         for (int j = 0; j < noccB; j++) {
              oo_idxAB->set(i, j, j + (i * noccB));
         }
    }

    // All oo
    oo_idxBB = SharedTensor2i(new Tensor2i("oo_idx", noccB, noccB));
    for (int i = 0; i < noccB; i++) {
         for (int j = 0; j < noccB; j++) {
              oo_idxBB->set(i, j, j + (i * noccB));
         }
    }

    // Active Ov
    ia_idxAB = SharedTensor2i(new Tensor2i("Ia_idx", naoccA, navirB));
    for (int i = 0; i < naoccA; i++) {
         for (int a = 0; a < navirB; a++) {
              ia_idxAB->set(i, a, a + (i * navirB));
         }
    }

    // Active oV
    ia_idxBA = SharedTensor2i(new Tensor2i("iA_idx", naoccB, navirA));
    for (int i = 0; i < naoccB; i++) {
         for (int a = 0; a < navirA; a++) {
              ia_idxBA->set(i, a, a + (i * navirA));
         }
    }

    // Active ov
    ia_idxBB = SharedTensor2i(new Tensor2i("ia_idx", naoccB, navirB));
    for (int i = 0; i < naoccB; i++) {
         for (int a = 0; a < navirB; a++) {
              ia_idxBB->set(i, a, a + (i * navirB));
         }
    }

    // All Ov
    ov_idxAB = SharedTensor2i(new Tensor2i("Ov_idx", noccA, nvirB));
    for (int i = 0; i < noccA; i++) {
         for (int a = 0; a < nvirB; a++) {
              ov_idxAB->set(i, a, a + (i * nvirB));
         }
    }

    // All ov
    ov_idxBB = SharedTensor2i(new Tensor2i("ov_idx", noccB, nvirB));
    for (int i = 0; i < noccB; i++) {
         for (int a = 0; a < nvirB; a++) {
              ov_idxBB->set(i, a, a + (i * nvirB));
         }
    }

    // Active Vo
    ai_idxAB = SharedTensor2i(new Tensor2i("Ai_idx", navirA, naoccB));
    for (int a = 0; a < navirA; a++) {
         for (int i = 0; i < naoccB; i++) {
              ai_idxAB->set(a, i, i + (a * naoccB));
         }
    }

    // Active vO
    ai_idxBA = SharedTensor2i(new Tensor2i("aI_idx", navirB, naoccA));
    for (int a = 0; a < navirB; a++) {
         for (int i = 0; i < naoccA; i++) {
              ai_idxBA->set(a, i, i + (a * naoccA));
         }
    }

    // Active vo
    ai_idxBB = SharedTensor2i(new Tensor2i("ai_idx", navirB, naoccB));
    for (int a = 0; a < navirB; a++) {
         for (int i = 0; i < naoccB; i++) {
              ai_idxBB->set(a, i, i + (a * naoccB));
         }
    }

    // All Vo
    vo_idxAB = SharedTensor2i(new Tensor2i("Vo_idx", nvirA, noccB));
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccB; i++) {
              vo_idxAB->set(a, i, i + (a * noccB));
         }
    }

    // All vo
    vo_idxBB = SharedTensor2i(new Tensor2i("vo_idx", nvirB, noccB));
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              vo_idxBB->set(a, i, i + (a * noccB));
         }
    }

    // Active Vv
    ab_idxAB = SharedTensor2i(new Tensor2i("Ab_idx", navirA, navirB));
    for (int a = 0; a < navirA; a++) {
         for (int b = 0; b < navirB; b++) {
              ab_idxAB->set(a, b, b + (a * navirB));
         }
    }

    // Active vV
    ab_idxBA = SharedTensor2i(new Tensor2i("aB_idx", navirB, navirA));
    for (int a = 0; a < navirB; a++) {
         for (int b = 0; b < navirA; b++) {
              ab_idxBA->set(a, b, b + (a * navirA));
         }
    }

    // Active vv
    ab_idxBB = SharedTensor2i(new Tensor2i("ab_idx", navirB, navirB));
    for (int a = 0; a < navirB; a++) {
         for (int b = 0; b < navirB; b++) {
              ab_idxBB->set(a, b, b + (a * navirB));
         }
    }

    // All Vv
    vv_idxAB = SharedTensor2i(new Tensor2i("Vv_idx", nvirA, nvirB));
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirB; b++) {
              vv_idxAB->set(a, b, b + (a * nvirB));
         }
    }

    // All vv
    vv_idxBB = SharedTensor2i(new Tensor2i("vv_idx", nvirB, nvirB));
    for (int a = 0; a < nvirB; a++) {
         for (int b = 0; b < nvirB; b++) {
              vv_idxBB->set(a, b, b + (a * nvirB));
         }
    }
 }// end if (reference_ == "UNRESTRICTED")
}//


//====================================
//   so_pair_idx (unpacked)
//====================================
int DFOCC::so_pair_idx(int i, int j)
{
   int value = j + (i * nso_);
   return value;
}

//====================================
//   mo_pair_idx (unpacked)
//====================================
int DFOCC::mo_pair_idx(int i, int j)
{
   int value = j + (i * nmo_);
   return value;
}

//====================================
//   oo_pair_idxAA (all OO)
//====================================
int DFOCC::oo_pair_idxAA(int i, int j)
{
   int value = j + (i * noccA);
   return value;
}

//====================================
//   ij_pair_idxAA (active OO)
//====================================
int DFOCC::ij_pair_idxAA(int i, int j)
{
   int value = j + (i * naoccA);
   return value;
}

//====================================
//   vv_pair_idxAA (all VV)
//====================================
int DFOCC::vv_pair_idxAA(int a, int b)
{
   int value = b + (a * nvirA);
   return value;
}

//====================================
//   ab_pair_idxAA (active VV)
//====================================
int DFOCC::ab_pair_idxAA(int a, int b)
{
   int value = b + (a * navirA);
   return value;
}

//====================================
//   ov_pair_idxAA (all OV)
//====================================
int DFOCC::ov_pair_idxAA(int i, int a)
{
   int value = a + (i * nvirA);
   return value;
}

//====================================
//   ia_pair_idxAA (active OV)
//====================================
int DFOCC::ia_pair_idxAA(int i, int a)
{
   int value = a + (i * navirA);
   return value;
}

//====================================
//   vo_pair_idxAA (all VVO
//====================================
int DFOCC::vo_pair_idxAA(int a, int i)
{
   int value = i + (a * noccA);
   return value;
}

//====================================
//   ai_pair_idxAA (active VVO
//====================================
int DFOCC::ai_pair_idxAA(int a, int i)
{
   int value = i + (a * naoccA);
   return value;
}

//====================================
//   MO Rotation Block
//====================================
int DFOCC::get_rotation_block(string rotblock)
{
   int index;

   if (rotblock == "VO")                         index = 1;
   else if (rotblock == "VO_AOCCFC")             index = 2;
   else if (rotblock == "VO_AOCCFC_FVAVIR")      index = 3;
   else                                          index = 4;     // NULL

   return(index);
}

}} // End Namespaces


