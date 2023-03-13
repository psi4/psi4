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

#ifndef PSI4_SRC_DLPNO_MP2_H_
#define PSI4_SRC_DLPNO_MP2_H_

#include "sparse.h"

#include "psi4/libmints/wavefunction.h"

#include <map>
#include <tuple>
#include <string>
#include <unordered_map>

namespace psi {
namespace dlpno {

enum VirtualStorage { DIRECT, DF_STORE, STORE };

// Equations refer to Pinski et al. (JCP 143, 034108, 2015; DOI: 10.1063/1.4926879)

class DLPNOMP2 : public Wavefunction {
   protected:

    /// threshold for PAO domain size
    double T_CUT_DO_;

    /// threshold for PNO truncation
    double T_CUT_PNO_;

    /// threshold for TNO truncation
    double T_CUT_TNO_;

    /// auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;
    SharedMatrix full_metric_;

    /// localized molecular orbitals (LMOs)
    std::map<std::string, SharedMatrix> lmo_matrices_; ///<A lookup table of each LMO matrix
    SharedMatrix C_lmo_;
    SharedMatrix F_lmo_;

    /// projected atomic orbitals (PAOs)
    std::map<std::string, SharedMatrix> pao_matrices_; ///<A lookup table of each PAO matrix
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

    /// LMO/PAO three-index integrals
    std::vector<SharedMatrix> qia_;

    /// pair natural orbitals (PNOs)
    std::map<std::string, std::vector<SharedMatrix>> pno_matrices_; ///<A lookup table of each PNO matrix
    std::vector<SharedMatrix> K_iajb_;  ///< exchange operators (i.e. (ia|jb) integrals)
    std::vector<SharedMatrix> T_iajb_;  ///< amplitudes
    std::vector<SharedMatrix> Tt_iajb_; ///< antisymmetrized amplitudes
    std::vector<SharedMatrix> T_ia_;   ///< singles amplitudes
    std::vector<SharedMatrix> X_pno_;   ///< global PAO -> canonical PNO transforms
    std::vector<SharedVector> e_pno_;   ///< PNO orbital energies
    std::vector<SharedMatrix> D_ij_; ///< pair densities
    std::vector<int> n_pno_;       ///< number of pnos
    std::vector<double> de_pno_;   ///< PNO truncation energy error
    std::vector<double> de_pno_os_;   ///< opposite-spin contributions to de_pno_
    std::vector<double> de_pno_ss_;   ///< same-spin contributions to de_pno_
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_mn_; ///< pno overlaps

    /// triplet natural orbitals (TNOs)
    std::vector<SharedMatrix> W_iajbkc_; ///< W3 intermediate for each lmo triplet
    std::vector<SharedMatrix> T_iajbkc_; ///< Triples amplitude for each lmo triplet
    std::vector<SharedMatrix> X_tno_; ///< global PAO -> canonical TNO transforms
    std::vector<SharedVector> e_tno_; ///< TNO orbital energies
    std::vector<SharedMatrix> denom_ijk_; /// (eps[a_ijk] + eps[b_ijk] + eps[c_ijk] - F[i,i] - F[j,j] - F[k,k])
    std::vector<int> n_tno_; ///<number of tnos per triplet domain
    std::vector<SharedMatrix> S_ijk_ii_; ///< tno/diagonal pno overlaps
    std::vector<SharedMatrix> S_ijk_ij_; ///< tno/pno overlaps
    std::vector<std::vector<SharedMatrix>> S_ijk_il_; ///< tno/pno overlaps
    std::vector<std::vector<SharedMatrix>> S_ijk_ljk_; ///< tno overlaps

    /// Integral Types
    std::vector<SharedMatrix> J_ijab_; /// (i j | a_ij b_ij)
    std::vector<SharedMatrix> K_mnij_; /// (m i | n j)
    std::vector<SharedMatrix> K_mbij_; /// (m i | b_ij j)
    std::vector<SharedMatrix> K_maef_; /// (m e_mm | a_mm f_mm)
    std::vector<std::vector<SharedMatrix>> Qab_ij_; // (q_ij | a_ij b_ij)
    std::vector<SharedMatrix> K_abef_; /// (a_ij e_ij | b_ij f_ij)
    std::vector<SharedMatrix> L_iajb_; /// 2.0 * (i a_ij | j b_ij) - (i b_ij | j a_ij)
    std::vector<SharedMatrix> Lt_iajb_; /// 2.0 * (i a_ij | j b_ij) - (i j | b_ij a_ij)

    // final energies
    double de_dipole_; ///< energy correction for distant (LMO, LMO) pairs
    double de_pno_total_; ///< energy correction for PNO truncation
    double de_pno_total_os_; ///< energy correction for PNO truncation
    double de_pno_total_ss_; ///< energy correction for PNO truncation
    double e_lmp2_; ///< raw (uncorrected) local MP2 correlation energy
    double e_lmp2_ss_; ///< same-spin component of e_lmp2_
    double e_lmp2_os_; ///< opposite-spin component of e_lmp2_
    double e_lccsd_; ///< raw (uncorrected) local CCSD correlation energy
    double e_lccsd_t_; ///< local (T) correlation energy

    // => Sparse Maps <= //
    std::map<std::string, SparseMap> sparse_maps_; ///< A lookup table of each SparseMap

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
    std::unordered_map<int, int> i_j_k_to_ijk_;  ///< LMO indices (i, j, k) to significant LMO triplet index (ijk)
    std::vector<std::tuple<int, int, int>> ijk_to_i_j_k_;  ///< LMO triplet index (ijk) to all three LMO indices (i, j, k)

    // LMO Pair Domains
    SparseMap lmopair_to_ribfs_; ///< which aux BFs are needed for density-fitting a pair of LMOs?
    SparseMap lmopair_to_riatoms_; ///< aux BFs on which atoms are needed for density-fitting a pair of LMOs?
    SparseMap lmopair_to_paos_; ///< which PAOs span the virtual space of a pair of LMOs?
    SparseMap lmopair_to_paoatoms_; ///< PAOs on which atoms span the virtual space of a pair of LMOs?
    SparseMap lmopair_to_lmos_; ///< Which LMOs "interact" with an LMO pair (determined by DOI integrals)
    // LMO triples domains
    SparseMap lmotriplet_to_paos_;
    // Useful for generating DF integrals
    std::vector<std::vector<std::vector<int>>> lmopair_lmo_to_riatom_lmo_;
    std::vector<std::vector<std::vector<int>>> lmopair_pao_to_riatom_pao_;

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

    /// CC Integrals
    VirtualStorage virtual_storage_;

    // LMO/LMO ERIs
    std::vector<SharedMatrix> qij_;
    // PAO/PAO ERIs
    std::vector<SharedMatrix> qab_;

    void common_init();

    std::pair<SharedMatrix, SharedVector> orthocanonicalizer(SharedMatrix S, SharedMatrix F);

    /// Form LMOs, PAOs, etc.
    void setup_orbitals();
    
    /// Compute differential overlap integrals between LMO/LMO and LMO/PAO pairs (EQ 4), DOI_ij and DOI_iu
    void compute_overlap_ints();

    /// Compute approximate MP2 pair energies for distant LMOs using dipole integrals (EQ 17)
    void compute_dipole_ints();

    /// Use dipole and overlap integrals to assess sparsity relationships between LMOs and estimate
    /// energy contribution of weakly interacting LMO pairs. Additionally, the overlap integrals
    /// and LMO sparsity are used to construct domains of PAOs, RI basis functions, and orbital
    /// basis functions for each LMO. These domains are necessary for efficient evaluation of
    /// three-index integrals.
    void prep_sparsity();

    /// Compute three-index integrals in LMO/PAO basis with linear scaling (EQ 11)
    void compute_df_ints();

    /// form pair exch operators (EQ 15) and SC amplitudes (EQ 18); transform to PNO basis
    void pno_transform();

    /// compute PNO/PNO overlap matrices
    void compute_pno_overlaps();

    /// compute MP2 correlation energy w/ current amplitudes (EQ 14)
    double compute_iteration_energy(const std::vector<SharedMatrix> &R_iajb);

    /// iteratively solve local MP2 equations  (EQ 13)
    void lmp2_iterations();

    // => CC intermediates <= //

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
    void print_aux_domains();
    void print_pao_domains();
    void print_lmo_domains();
    void print_aux_pair_domains();
    void print_pao_pair_domains();
    void print_integral_sparsity();
    void print_results();
    
    void print_ccsd_results();
    void print_ccsd_t_results();

    /// Create TNOs (Triples Natural Orbitals) for DLPNO-CCSD(T)
    void tno_transform();
    /// Compute TNO/TNO overlap matrices
    void compute_tno_overlaps();

    /// compute W3 intermediate (for DLPNO-CCSD(T))
    void compute_W_iajbkc();
    /// compute (T) energy
    double compute_t_energy();
    /// L_CCSD(T0) energy
    void compute_lccsd_t0();
    /// L_CCSD(T) iterations
    void lccsd_t_iterations();

    /// Create lookup tables for all sparsity information (LMO, PAO, PNOs, SparseMaps)
    void store_information();

    /// A function to estimate integral memory costs
    void estimate_memory();

    /// A function to compute Qij integrals
    void compute_qij();
    /// A function to compute Qab integrals
    void compute_qab();

    /// compute CC PNO overlaps
    void compute_cc_pno_overlaps();

    /// Returns the appropriate overlap matrix given two LMO pairs
    inline SharedMatrix S_PNO(int ij, int mn);

    /// Compute four-center integrals for CC computations
    void compute_cc_ints();

   public:
    DLPNOMP2(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOMP2() override;

    /// Gets LMO matrix from LMO lookup table
    SharedMatrix get_lmo_matrix(std::string key);
    /// Gets PAO matrix from PAO lookup table
    SharedMatrix get_pao_matrix(std::string key);
    /// Gets PNO matrix from PNO lookup table
    std::vector<SharedMatrix> get_pno_matrix(std::string key);
    /// Gets specific sparse map by key
    SparseMap get_sparse_map(std::string key);
    /// Returns PNO virtual energies per pair ij
    std::vector<SharedVector> eps_pno() { return e_pno_; }
    /// Returns X_tno matrix
    std::vector<SharedMatrix> get_X_tno() { return X_tno_; }
    /// Returns TNO virtual energies for triplet ijk
    std::vector<SharedVector> eps_tno() { return e_tno_; }

    /// LMO indices (i, j) to significant LMO pair index (ij); insignificant (i, j) maps to -1
    std::vector<std::vector<int>> i_j_to_ij() { return i_j_to_ij_; }
    /// LMO pair index (ij) to both LMO indices (i, j)
    std::vector<std::pair<int,int>> ij_to_i_j() { return ij_to_i_j_; }
    /// LMO pair index (ij) to LMO pair index (ji)
    std::vector<int> ij_to_ji() { return ij_to_ji_; }

    /// LMO indices (i, j, k) to significant LMO triplet index (ijk)
    std::unordered_map<int, int> i_j_k_to_ijk() { return i_j_k_to_ijk_; }
    /// LMO triplet index (ijk) to all three LMO indices (i, j, k)
    std::vector<std::tuple<int, int, int>> ijk_to_i_j_k() { return ijk_to_i_j_k_; }

    /// Geta qia matrix
    std::vector<SharedMatrix> get_qia();
    /// Gets qij matrix
    std::vector<SharedMatrix> get_qij();
    /// Gets qab matrix
    std::vector<SharedMatrix> get_qab();

    /// Returns LMP2 Correlation Energy
    double e_lmp2() { return e_lmp2_; }
    /// Returns LMO Truncation Correction
    double de_lmo_total() { return de_dipole_; }
    /// Returns PNO Truncation Correction
    double de_pno_total() { return de_pno_total_; }

    /// Computes information for DLPNO-(T)
    void compute_triples_info();

    double compute_energy() override;

};

}  // namespace dlpno
}  // namespace psi

#endif //PSI4_SRC_DLPNO_MP2_H_
