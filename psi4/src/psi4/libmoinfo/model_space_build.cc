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

#include "model_space.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdio>
namespace psi {

void ModelSpace::build() {
    /********************************************************
      Generate all the model space Slater determinants
      - docc MOs are doubly    occupied
      - actv MOs are partially occupied
      - the generalized occupied orbital indexing (docc + actv)
        is assumed (see the MOInfo class)
      - external orbitals are omitted
      - the table below illustrates the orbital ordering

      ----------------------------
      |0 1 2 3 4 5| 6 7 8 9 10 11|
      ----------------------------
      |  Irrep 0  |   Irrep 1    |
      | docc |act | docc   | act |
      |     occ   |       occ    |
      |--------------------------|

     ********************************************************/

    int nirreps = moinfo_obj->get_nirreps();
    int nactv = moinfo_obj->get_nactv();
    int nactive_ael = moinfo_obj->get_nactive_ael();
    int nactive_bel = moinfo_obj->get_nactive_bel();

    std::vector<int> docc = moinfo_obj->get_docc();
    std::vector<int> actv = moinfo_obj->get_actv();

    std::vector<bool> docc_bits;

    //  Set up the doubly occupied part of all the
    //  determinants in the model space
    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < docc[h]; ++i) docc_bits.push_back(true);
        for (int i = 0; i < actv[h]; ++i) docc_bits.push_back(false);
    }

    //  Set up the a vectors containing the active orbitals
    //  in and their symmetry
    std::vector<int> active, active_sym;
    int mo_index = 0;
    for (int h = 0; h < nirreps; ++h) {
        mo_index += docc[h];
        for (int i = 0; i < actv[h]; ++i) {
            active.push_back(mo_index);
            active_sym.push_back(h);
            mo_index++;
        }
    }

    // Generate combinations of alpha and beta orbitals
    std::vector<std::vector<int> > alfa_combinations, beta_combinations;
    generate_combinations(nactv, nactive_ael, alfa_combinations);
    generate_combinations(nactv, nactive_bel, beta_combinations);

    if (alfa_combinations.size() == 0) alfa_combinations.push_back(std::vector<int>(0));
    if (beta_combinations.size() == 0) beta_combinations.push_back(std::vector<int>(0));

    for (size_t a = 0; a < alfa_combinations.size(); ++a) {
        for (size_t b = 0; b < beta_combinations.size(); ++b) {
            int alfa_sym = 0;  // Symmetry of the alfa string
            int beta_sym = 0;  // Symmetry of the beta string

            // Fill the doubly occupied MOs
            std::vector<bool> alfa_bits = docc_bits;
            std::vector<bool> beta_bits = docc_bits;

            // Fill the alpha active orbitals
            for (int i = 0; i < nactive_ael; ++i) {
                alfa_bits[active[alfa_combinations[a][i]]] = true;
                alfa_sym = alfa_sym ^ active_sym[alfa_combinations[a][i]];
            }

            // Fill the beta active orbitals
            for (int i = 0; i < nactive_bel; ++i) {
                beta_bits[active[beta_combinations[b][i]]] = true;
                beta_sym = beta_sym ^ active_sym[beta_combinations[b][i]];
            }

            // Test the symmetry of the determinant generated
            int sym = alfa_sym ^ beta_sym;
            if (sym == wfn_sym) {
                SlaterDeterminant det(alfa_sym, beta_sym, alfa_bits, beta_bits);
                determinants.push_back(det);
            }
        }
    }

    if (determinants.size() == 0) {
        std::string message =
            "\n\n  No reference found in the model space"
            "\n  Please check the following:"
            "\n  1) Definition of FROZEN_DOCC, RESTRICTED_DOCC, ACTIVE, and FROZEN_UOCC"
            "\n  2) Symmetry of the wavefunction"
            "\n  3) Charge and multiplicity"
            "\n\n  Ending the computation.\n";
        outfile->Printf(message);

        throw std::logic_error(message);
    }
}

void ModelSpace::classify() {
    for (size_t mu = 0; mu < determinants.size(); ++mu) {
        if (determinants[mu].is_closed_shell()) {
            closed_to_all.push_back(mu);
        } else {
            opensh_to_all.push_back(mu);
        }
    }
}

}  // namespace psi

//  // Closed-shell determinant
//  if(alfa_bits == beta_bits){
//    closed_to_all.push_back(mu);
//    unique_to_all.push_back(mu);
//  }else{
//    // Open-shell determinant
//    bool add_it = true;
//    int  spin_mirror = mu;
//    if(true /*options_get_bool("USE_SPIN_SYMMETRY")*/){
//      // Check if this is a spin-flipped determinant
//      for(int mu = 0; mu < reference; ++det){
//        if(references[ref].is_spin_flipped(det)){
//          add_it      = false;
//          spin_mirror = ref;
//        }
//      }
//      if(add_it){
//        unique_open_shell_refs.push_back(references.size());
//        unique_refs.push_back(references.size());
//        all_refs.push_back(references.size());
//      }else{
//        all_refs.push_back(spin_mirror);
//      }
//    }
//  }
//}
//}
//  }else{
//    /********************************************************
//      Set up the doubly occupied part of all the
//      determinants in the model space
//    ********************************************************/
//    index = 0;
//    for(int h=0;h<nirreps;h++){
//      for(int i=0;i<docc[h] + actv_docc[h];i++){
//        docc_det.set(index);
//        docc_det.set(index + nall);
//        index++;
//      }
//      index += actv[h] - actv_docc[h];
//      index += extr[h];
//    }
//    closed_shell_refs.push_back(references.size());
//    unique_refs.push_back(references.size());
//    all_refs.push_back(references.size());
//    references.push_back(docc_det);
//  }
