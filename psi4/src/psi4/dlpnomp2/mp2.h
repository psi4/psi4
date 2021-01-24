/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace dlpnomp2 {

// References to DiStasio are to J Comput Chem 28: 839–856, 2007; DOI: 10.1002/jcc.20604

class DLPNOMP2 : public Wavefunction {
   protected:
    double de_dipole_;
    double de_pno_total_;
   
    // Auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;

    SharedMatrix C_lmo;
    SharedMatrix F_lmo;
    SharedMatrix S_lmo;

    SharedMatrix C_pao;
    SharedMatrix F_pao;
    SharedMatrix S_pao;

    SharedMatrix DOI_ij;
    SharedMatrix DOI_iu;

    SharedMatrix e_actual;
    SharedMatrix e_linear;
    SharedMatrix lmo_distances;

    std::vector<SharedMatrix> qia;

    SharedMatrix full_metric;

    std::vector<std::vector<int>> atom_to_bf_;
    std::vector<int> bf_to_atom_;

    std::vector<std::vector<int>> atom_to_ribf_;
    std::vector<int> ribf_to_atom_;

    std::vector<std::vector<int>> atom_to_shell_;
    std::vector<std::vector<int>> atom_to_rishell_;

    // aux domains
    // sparse map from local MOs to auxiliary basis functions in the LMO's domain
    std::vector<std::vector<int>> lmo_to_ribfs_;
    // same  as above, but map to atom indices instead of bf indices
    std::vector<std::vector<int>> lmo_to_riatoms_;

    // pao domains
    // sparse map from LMO [i] to list of PAOs [u] in the domain of [i].
    std::vector<std::vector<int>> lmo_to_paos_;
    // same  as above, but map to atom indices instead of pao indices
    std::vector<std::vector<int>> lmo_to_paoatoms_;

    // lmo domains
    // map from LMO indices (i, j) to LMO pair index (ij)
    std::vector<std::vector<int>> i_j_to_ij;
    // map from LMO pair index (ij) to both LMO indices (i, j)
    std::vector<std::pair<int,int>> ij_to_i_j;
    // map from LMO pair index (ij) to LMO pair index (ji)
    std::vector<int> ij_to_ji;



    // sparse map from LMO pair ij to all PAO's u in the pair domain of i and j.
    std::vector<std::vector<int>> lmopair_to_paos;
    std::vector<std::vector<int>> lmopair_to_paoatoms;

    // sparse map from LMO pair ij to all AUXBF u (or atom a) in the pair domain of i and j.
    std::vector<std::vector<int>> lmopair_to_ribfs;
    std::vector<std::vector<int>> lmopair_to_riatoms;



    // extended fitting map (include fitting domains of all local MOs in your domain)
    std::vector<std::vector<int>> lmo_to_riatoms_ext;

    // target of the integral transformation
    std::vector<std::vector<int>> riatom_to_lmos_ext;
    std::vector<std::vector<int>> riatom_to_paos_ext;

    // We'll use these maps to screen the local MO transform: (mn|Q) * C_mi -> (in|Q)
    std::vector<std::vector<int>> riatom_to_atoms1;
    std::vector<std::vector<int>> riatom_to_shells1;
    std::vector<std::vector<int>> riatom_to_bfs1;

    // We'll use these maps to screen the projected AO transform: (mn|Q) * C_nu -> (mu|Q) 
    std::vector<std::vector<int>> riatom_to_atoms2;
    std::vector<std::vector<int>> riatom_to_shells2;
    std::vector<std::vector<int>> riatom_to_bfs2;

    // Need dense versions of previous maps for quick lookup
    // arr[riatom][lmo] is the index of lmo in riatom_to_lmos_ext[riatom] (if present), else -1
    std::vector<std::vector<int>> riatom_to_lmos_ext_dense;
    std::vector<std::vector<bool>> riatom_to_atoms1_dense;
    std::vector<std::vector<bool>> riatom_to_atoms2_dense;


    std::vector<SharedMatrix> K_iajb;  // exchange operators (i.e. (ia|jb) integrals)
    std::vector<SharedMatrix> T_iajb;  // amplitudes
    std::vector<SharedMatrix> Tt_iajb; // antisymmetrized amplitudes
    std::vector<SharedMatrix> X_pno;   // global PAOs -> canonical PNOs
    std::vector<SharedVector> e_pno;   // PNO orbital energies

    std::vector<int> n_pno;         // number of pnos
    std::vector<double> de_pno;   // PNO truncation error

    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj; // ij kj
    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik; // ij ik

    std::vector<double> iter_energies;
    std::vector<double> iter_delta_energies;
    std::vector<double> iter_max_residuals;

    void common_init();

    void overlap_ints();
    void dipole_ints();
    void sparsity_prep();
    void df_ints();
    void pno_transform();
    void pno_overlaps();
    double eval_amplitudes();
    void lmp2_iterations();


    void print_header();

    void print_aux_domains();
    void print_pao_domains();
    void print_lmo_domains();

    void print_aux_pair_domains();
    void print_pao_pair_domains();

    void print_integral_sparsity();

    void print_results();


   public:
    DLPNOMP2(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOMP2() override;

    double compute_energy();
};

}  // namespace dlpnomp2
}  // namespace psi
