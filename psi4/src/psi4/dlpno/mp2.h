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

#ifndef PSI4_SRC_DLPNO_MP2_H_
#define PSI4_SRC_DLPNO_MP2_H_

#include "sparse.h"

#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace dlpno {

// Equations refer to Pinski et al. (JCP 143, 034108, 2015; DOI: 10.1063/1.4926879)

class DLPNOMP2 : public Wavefunction {
   protected:

    /// threshold for PAO domain size
    double T_CUT_DO;

    /// threshold for PNO truncation
    double T_CUT_PNO;

    /// auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;
    SharedMatrix full_metric;

    /// localized molecular orbitals (LMOs)
    SharedMatrix C_lmo;
    SharedMatrix F_lmo;

    /// projected atomic orbitals (PAOs)
    SharedMatrix C_pao;
    SharedMatrix F_pao;
    SharedMatrix S_pao;

    /// differential overlap integrals (EQ 4)
    SharedMatrix DOI_ij; // LMO/LMO
    SharedMatrix DOI_iu; // LMO/PAO

    // approximate LMO/LMO pair energies from dipole integrals (EQ 17)
    SharedMatrix e_actual; ///< actual approximate pair energy
    SharedMatrix e_linear; ///< upper bound to approximate pair energies

    /// LMO/PAO three-index integrals
    std::vector<SharedMatrix> qia;

    /// pair natural orbitals (PNOs)
    std::vector<SharedMatrix> K_iajb;  ///< exchange operators (i.e. (ia|jb) integrals)
    std::vector<SharedMatrix> T_iajb;  ///< amplitudes
    std::vector<SharedMatrix> Tt_iajb; ///< antisymmetrized amplitudes
    std::vector<SharedMatrix> R_iajb;  ///< LMP2 residuals
    std::vector<SharedMatrix> X_pno;   ///< global PAO -> canonical PNO transforms
    std::vector<SharedVector> e_pno;   ///< PNO orbital energies
    std::vector<int> n_pno;       ///< number of pnos
    std::vector<double> de_pno;   ///< PNO truncation energy error
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik; ///< pnooverlaps

    // final energies
    double de_dipole_; ///< energy correction for distant (LMO, LMO) pairs
    double de_pno_total_; ///< energy correction for PNO truncation
    double e_lmp2_; ///< raw (uncorrected) local MP2 correlation energy

    // => Sparse Maps <= //

    // orbital / aux bases
    SparseMap atom_to_bf_; ///< which orbital BFs are on a given atom?
    SparseMap atom_to_ribf_; ///< which aux BFs are on a given atom?
    SparseMap atom_to_shell_; ///< which orbital basis shells are on a given atom?
    SparseMap atom_to_rishell_; ///< which aux basis shells are on a given atom?

    // LMO domains
    SparseMap lmo_to_ribfs_; ///< which aux BFs are needed for density-fitting a LMO?
    SparseMap lmo_to_riatoms_; ///< aux BFs on which atoms are needed for density-fitting a LMO?
    SparseMap lmo_to_paos_; ///< which PAOs span the virtual space of a LMO?
    SparseMap lmo_to_paoatoms_; ///< PAOs on which atoms span the virtual space of a LMO?
    std::vector<std::vector<int>> i_j_to_ij; ///< LMO indices (i, j) to LMO pair index (ij)
    std::vector<std::pair<int,int>> ij_to_i_j; ///< LMO pair index (ij) to both LMO indices (i, j)
    std::vector<int> ij_to_ji; ///< LMO pair index (ij) to LMO pair index (ji)

    // LMO Pair Domains
    SparseMap lmopair_to_ribfs; ///< which aux BFs are needed for density-fitting a pair of LMOs?
    SparseMap lmopair_to_riatoms; ///< aux BFs on which atoms are needed for density-fitting a pair of LMOs?
    SparseMap lmopair_to_paos; ///< which PAOs span the virtual space of a pair of LMOs?
    SparseMap lmopair_to_paoatoms; ///< PAOs on which atoms span the virtual space of a pair of LMOs?

    // Extended LMO Domains 
    SparseMap lmo_to_riatoms_ext; ///< aux BFs on which atoms are needed for density-fitting a LMO and all connected LMOs
    SparseMap riatom_to_lmos_ext; ///< the extended DF domains of which LMOs include aux BFs on an atom
    SparseMap riatom_to_paos_ext; ///< the extended DF domains of which PAOs include aux BFs on an atom
    SparseMap riatom_to_atoms1; ///< orbital BFs on which atoms are needed for DF int transform (first index)
    SparseMap riatom_to_shells1; ///< which shells of orbital BFs are needed for DF int transform (first index)
    SparseMap riatom_to_bfs1; ///< which orbital BFs are needed for DF int transform (first index)
    SparseMap riatom_to_atoms2; ///< orbital BFs on which atoms are needed for DF int transform (second index)
    SparseMap riatom_to_shells2; ///< which shells of orbital BFs are needed for DF int transform (second index)
    SparseMap riatom_to_bfs2; ///< which orbital BFs are needed for DF int transform (second index)

    // Dense analogues of some sparse maps for quick lookup
    std::vector<std::vector<int>> riatom_to_lmos_ext_dense;
    std::vector<std::vector<bool>> riatom_to_atoms1_dense;
    std::vector<std::vector<bool>> riatom_to_atoms2_dense;

    void common_init();

    std::pair<SharedMatrix, SharedVector> orthocanonicalizer(SharedMatrix S, SharedMatrix F);

    /// form LMOs, PAOs, etc.
    void setup();
    
    /// compute differential overlap integrals between LMO/LMO and LMO/PAO pairs (EQ 4)
    void overlap_ints();

    /// compute approximate MP2 pair energies for distant LMOs using dipole integrals (EQ 17)
    void dipole_ints();

    /// assess sparsity relationships between LMOs, PAOs, RI basis, and orbital basis
    void sparsity_prep();

    /// compute three-index integrals in LMO/PAO basis with linear scaling (EQ 11)
    void df_ints();

    /// form pair exch operators (EQ 15) and SC amplitudes (EQ 18); transform to PNO basis
    void pno_transform();

    /// compute PNO/PNO overlap matrices
    void pno_overlaps();

    /// compute MP2 correlation energy w/ current amplitudes (EQ 14)
    double eval_amplitudes();

    /// iteratively solve local MP2 equations  (EQ 13)
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

}  // namespace dlpno
}  // namespace psi

#endif //PSI4_SRC_DLPNO_MP2_H_
