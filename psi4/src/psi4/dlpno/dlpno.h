/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#ifndef PSI4_SRC_DLPNO_H_
#define PSI4_SRC_DLPNO_H_

#include "sparse.h"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"

#include <map>
#include <tuple>
#include <string>
#include <unordered_map>

namespace psi {
namespace dlpno {

enum VirtualStorage { CORE, DIRECT };

enum AlgorithmType { MP2, CCSD };

// Equations refer to Pinski et al. (JCP 143, 034108, 2015; DOI: 10.1063/1.4926879)

class DLPNOBase : public Wavefunction {
    protected:
      /// what quantum chemistry module are we running
      AlgorithmType algorithm_;

      /// threshold for PAO domain size
      double T_CUT_DO_;

      /// threshold for PNO truncation
      double T_CUT_PNO_;

      /// tolerance to separate pairs into CCSD and MP2 pairs
      double T_CUT_PAIRS_;

      /// tolerance for local density fitting (by Mulliken population)
      double T_CUT_MKN_;

      /// T_CUT_PNO scaling factor for diagonal PNOs
      double T_CUT_PNO_DIAG_SCALE_;

      

      /// auxiliary basis
      std::shared_ptr<BasisSet> ribasis_;
      SharedMatrix full_metric_;

      /// localized molecular orbitals (LMOs)
      SharedMatrix C_lmo_;
      SharedMatrix F_lmo_;

      /// projected atomic orbitals (PAOs)
      SharedMatrix C_pao_;
      SharedMatrix F_pao_;
      SharedMatrix S_pao_;

      /// differential overlap integrals (EQ 4)
      SharedMatrix DOI_ij_; // LMO/LMO
      SharedMatrix DOI_iu_; // LMO/PAO

      // approximate LMO/LMO pair energies from dipole integrals (EQ 17)
      // used to screen out and estimate weakly interacting LMO/LMO pairs
      SharedMatrix dipole_pair_e_; ///< actual approximate pair energy (used in final energy calculation)
      SharedMatrix dipole_pair_e_bound_; ///< upper bound to approximate pair energies (used for screening)

      /// LMO/LMO three-index integrals
      std::vector<SharedMatrix> qij_;
      /// LMO/PAO three-index integrals
      std::vector<SharedMatrix> qia_;
      /// PAO/PAO three-index integrals
      std::vector<SharedMatrix> qab_;

      /// pair natural orbitals (PNOs)
      std::vector<SharedMatrix> K_iajb_;  ///< exchange operators (i.e. (ia|jb) integrals)
      std::vector<SharedMatrix> T_iajb_;  ///< amplitudes
      std::vector<SharedMatrix> Tt_iajb_; ///< antisymmetrized amplitudes
      std::vector<SharedMatrix> X_pno_;   ///< global PAO -> canonical PNO transforms
      std::vector<SharedVector> e_pno_;   ///< PNO orbital energies
      std::vector<int> n_pno_;       ///< number of pnos
      std::vector<double> de_pno_;   ///< PNO truncation energy error
      std::vector<double> de_pno_os_;   ///< opposite-spin contributions to de_pno_
      std::vector<double> de_pno_ss_;   ///< same-spin contributions to de_pno_

      /// pre-screening energies
      double de_dipole_; ///< energy correction for distant (LMO, LMO) pairs
      double de_lmp2_; ///< SC-LMP2 correction for weak pairs (only for CC)
      double de_pno_total_; ///< energy correction for PNO truncation
      double de_pno_total_os_; ///< energy correction for PNO truncation
      double de_pno_total_ss_; ///< energy correction for PNO truncation

      // => Sparse Maps <= //

      // orbital / aux bases
      SparseMap atom_to_bf_; ///< which orbital BFs are on a given atom?
      SparseMap atom_to_ribf_; ///< which aux BFs are on a given atom?
      SparseMap atom_to_shell_; ///< which orbital basis shells are on a given atom?
      SparseMap atom_to_rishell_; ///< which aux basis shells are on a given atom?

      // AO to LMO/PAO
      SparseMap lmo_to_bfs_;
      SparseMap lmo_to_atoms_;
      SparseMap pao_to_bfs_;
      SparseMap pao_to_atoms_;

      // LMO domains
      SparseMap lmo_to_ribfs_; ///< which aux BFs are needed for density-fitting a LMO?
      SparseMap lmo_to_riatoms_; ///< aux BFs on which atoms are needed for density-fitting a LMO?
      SparseMap lmo_to_paos_; ///< which PAOs span the virtual space of a LMO?
      SparseMap lmo_to_paoatoms_; ///< PAOs on which atoms span the virtual space of a LMO?
      std::vector<std::vector<int>> i_j_to_ij_; ///< LMO indices (i, j) to significant LMO pair index (ij); insignificant (i, j) maps to -1
      std::vector<std::pair<int,int>> ij_to_i_j_; ///< LMO pair index (ij) to both LMO indices (i, j)
      std::vector<int> ij_to_ji_; ///< LMO pair index (ij) to LMO pair index (ji)

      // LMO Pair Domains
      SparseMap lmopair_to_ribfs_; ///< which aux BFs are needed for density-fitting a pair of LMOs?
      SparseMap lmopair_to_riatoms_; ///< aux BFs on which atoms are needed for density-fitting a pair of LMOs?
      SparseMap lmopair_to_paos_; ///< which PAOs span the virtual space of a pair of LMOs?
      SparseMap lmopair_to_paoatoms_; ///< PAOs on which atoms span the virtual space of a pair of LMOs?
      SparseMap lmopair_to_lmos_; ///< Which LMOs "interact" with an LMO pair (determined by DOI integrals)

      // Extended LMO Domains 
      SparseMap lmo_to_riatoms_ext_; ///< aux BFs on which atoms are needed for density-fitting a LMO and all connected LMOs
      SparseMap riatom_to_lmos_ext_; ///< the extended DF domains of which LMOs include aux BFs on an atom
      SparseMap riatom_to_paos_ext_; ///< the extended DF domains of which PAOs include aux BFs on an atom
      SparseMap riatom_to_atoms1_; ///< orbital BFs on which atoms are needed for DF int transform (first index)
      SparseMap riatom_to_shells1_; ///< which shells of orbital BFs are needed for DF int transform (first index)
      SparseMap riatom_to_bfs1_; ///< which orbital BFs are needed for DF int transform (first index)
      SparseMap riatom_to_atoms2_; ///< orbital BFs on which atoms are needed for DF int transform (second index)
      SparseMap riatom_to_shells2_; ///< which shells of orbital BFs are needed for DF int transform (second index)
      SparseMap riatom_to_bfs2_; ///< which orbital BFs are needed for DF int transform (second index)

      // Dense analogues of some sparse maps for quick lookup
      std::vector<std::vector<int>> riatom_to_lmos_ext_dense_;
      std::vector<std::vector<int>> riatom_to_paos_ext_dense_;
      std::vector<std::vector<bool>> riatom_to_atoms1_dense_;
      std::vector<std::vector<bool>> riatom_to_atoms2_dense_;
      std::vector<std::vector<int>> lmopair_to_lmos_dense_;

      /// Useful for generating DF integrals
      std::vector<std::vector<std::vector<int>>> lmopair_lmo_to_riatom_lmo_;
      std::vector<std::vector<std::vector<int>>> lmopair_pao_to_riatom_pao_;

      void common_init();

      // Helper functions
      void C_DGESV_wrapper(SharedMatrix A, SharedMatrix B);

      std::pair<SharedMatrix, SharedVector> canonicalizer(SharedMatrix C, SharedMatrix F);
      std::pair<SharedMatrix, SharedVector> orthocanonicalizer(SharedMatrix S, SharedMatrix F);

      SharedVector flatten_mats(const std::vector<SharedMatrix>& mat_list);

      void copy_flat_mats(SharedVector flat, std::vector<SharedMatrix>& mat_list);

      /// Form LMOs, PAOs, etc.
      void setup_orbitals();

      /// Compute approximate MP2 pair energies for distant LMOs using dipole integrals (EQ 17)
      void compute_dipole_ints();

      /// Compute differential overlap integrals between LMO/LMO and LMO/PAO pairs (EQ 4), DOI_ij and DOI_iu
      void compute_overlap_ints();

      /// Use dipole and overlap integrals to assess sparsity relationships between LMOs and estimate
      /// energy contribution of weakly interacting LMO pairs. Additionally, the overlap integrals
      /// and LMO sparsity are used to construct domains of PAOs, RI basis functions, and orbital
      /// basis functions for each LMO. These domains are necessary for efficient evaluation of
      /// three-index integrals.
      void prep_sparsity();

      /// Compute the auxiliary metric (P|Q)
      void compute_metric();
      /// Compute three-index integrals in LMO/LMO basis with linear scaling
      void compute_qij();
      /// Compute three-index integrals in LMO/PAO basis with linear scaling (EQ 11)
      void compute_qia();
      /// Compute three-index integrals in PAO/PAO basis with linear scaling
      void compute_qab();

      /// form pair exch operators (EQ 15) and SC amplitudes (EQ 18); transform to PNO basis
      void pno_transform();

      /// Printing
      void print_aux_domains();
      void print_pao_domains();
      void print_lmo_domains();
      void print_aux_pair_domains();
      void print_pao_pair_domains();

    public:
      DLPNOBase(SharedWavefunction ref_wfn, Options& options);
      ~DLPNOBase() override;

      double compute_energy() override;
};

class DLPNOMP2 : public DLPNOBase {
   protected:
    /// PNO overlap integrals
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik_; ///< pno overlaps

    double e_lmp2_; ///< raw (uncorrected) local MP2 correlation energy
    double e_lmp2_ss_; ///< same-spin component of e_lmp2_
    double e_lmp2_os_; ///< opposite-spin component of e_lmp2_

    /// compute PNO/PNO overlap matrices for DLPNO-MP2
    void compute_pno_overlaps();
    
    /// compute MP2 correlation energy w/ current amplitudes (EQ 14)
    double compute_iteration_energy(const std::vector<SharedMatrix> &R_iajb);

    /// iteratively solve local MP2 equations  (EQ 13)
    void lmp2_iterations();

    void print_header();
    void print_results();
    void print_integral_sparsity();

   public:
    DLPNOMP2(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOMP2() override;

    double compute_energy() override;
};

class DLPNOCCSD : public DLPNOBase {
   protected:
    /// How to compute CC integrals
    VirtualStorage virtual_storage_;

    /// How much memory is used by storing each of the DF integral types
    size_t qij_memory_;
    size_t qia_memory_;
    size_t qab_memory_;

    /// PNO overlap integrals
    std::vector<std::vector<SharedMatrix>> S_pno_ij_mn_; ///< pno overlaps

    /// Coupled-cluster amplitudes
    std::vector<SharedMatrix> T_ia_; ///< singles amplitudes

    // => Strong and Weak Pair Info <=//

    std::vector<std::vector<int>> i_j_to_ij_strong_;
    std::vector<std::pair<int,int>> ij_to_i_j_strong_;
    std::vector<int> ij_to_ji_strong_;

    std::vector<std::vector<int>> i_j_to_ij_weak_;
    std::vector<std::pair<int,int>> ij_to_i_j_weak_;
    std::vector<int> ij_to_ji_weak_;

    // => CCSD Integrals <= //

    /// (4 occupied, 0 virtual)
    std::vector<SharedMatrix> K_mnij_; /// (m i | n j)
    /// (3 occupied, 1 virtual)
    std::vector<SharedMatrix> K_bar_; /// (m i | b_ij j) [aka K_bar]
    std::vector<SharedMatrix> L_bar_; /// 2.0 * K_mbij - K_mbji
    /// (2 occupied, 2 virtual)
    std::vector<SharedMatrix> J_ijab_; /// (i j | a_ij b_ij)
    std::vector<SharedMatrix> L_iajb_; /// 2.0 * (i a_ij | j b_ij) - (i b_ij | j a_ij)
    std::vector<SharedMatrix> M_iajb_; /// 2.0 * (i a_ij | j b_ij) - (i j | b_ij a_ij)
    /// (1 occupied, 3 virtual)
    std::vector<SharedMatrix> K_tilde_chem_; /// (i e_ij | a_ij f_ij) [aka K_tilde] (stored as (e, a*f)) [Chemist's Notation]
    std::vector<SharedMatrix> K_tilde_phys_; /// (i e_ij | a_ij f_ij) [aka K_tilde] (stored as (a, e*f)) [Physicist's Notation]
    std::vector<SharedMatrix> L_tilde_; /// 2.0 * K_tilde_chem - K_tilde_phys
    /// (0 occupied, 4 virtual)
    std::vector<std::vector<SharedMatrix>> Qab_ij_; // (q_ij | a_ij b_ij)

    double e_lccsd_; ///< raw (uncorrected) local CCSD correlation energy

    /// Returns the appropriate overlap matrix given two LMO pairs
    inline SharedMatrix S_PNO(const int ij, const int mn);

    /// Determine which pairs are strong and weak pairs
    void determine_strong_and_weak_pairs();

    /// compute PNO/PNO overlap matrices for DLPNO-CCSD
    void compute_pno_overlaps();

    // => Computing integrals <= //

    /// A function to estimate integral memory costs
    void estimate_memory();
    /// Compute four-center integrals for CC computations
    void compute_cc_integrals();

    // => CCSD intermediates <= //

    /// compute Fmi intermediate (Madriaga Eq. 40)
    SharedMatrix compute_Fmi(const std::vector<SharedMatrix>& tau_tilde);
    /// compute Fbe intermediate (of diagonal LMO pair ii) (Madriaga Eq. 39)
    std::vector<SharedMatrix> compute_Fbe(const std::vector<SharedMatrix>& tau_tilde);
    /// compute Fme intermediate (of diagonal LMO pair mm) (Madriaga Eq. 41)
    std::vector<SharedMatrix> compute_Fme();
    /// compute Wmnij intermediate (Madriaga Eq. 43)
    std::vector<SharedMatrix> compute_Wmnij(const std::vector<SharedMatrix>& tau);
    /// compute Wmbej intermediate (Madriaga Eq. 44)
    std::vector<SharedMatrix> compute_Wmbej(const std::vector<SharedMatrix>& tau_bar);
    /// compute Wmbje intermediate (Madriaga Eq. 45)
    std::vector<SharedMatrix> compute_Wmbje(const std::vector<SharedMatrix>& tau_bar);

    /// iteratively solve local CCSD equations
    void lccsd_iterations();

    void print_header();
    void print_results();
    void print_integral_sparsity();
    
   public:
    DLPNOCCSD(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSD() override;

    double compute_energy() override;
};

}  // namespace dlpno
}  // namespace psi

#endif //PSI4_SRC_DLPNO_MP2_H_
