/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"

#include <map>
#include <tuple>
#include <string>
#include <unordered_map>

#ifdef USING_Einsums
#include "Einsums/Tensor.hpp"
#include "Einsums/TensorAlgebra.hpp"
#include "Einsums/LinearAlgebra.hpp"
#include "Einsums/Profile.hpp"
#include "Einsums/TensorUtilities/RMSD.hpp"
#endif

namespace psi {
namespace dlpno {

enum class DLPNOMethod { MP2, CCSD, CCSD_T, CCSDT, CCSDT_Q, CCSDTQ };

/// Flatten an occupied-orbital triplet into a dense, size-safe lookup key.
/// Shared by the triples and quadruples implementations so unity builds do not
/// define two same-named anonymous-namespace helpers in one translation unit.
inline size_t triplet_key(int i, int j, int k, size_t nocc) {
    return (static_cast<size_t>(i) * nocc + j) * nocc + k;
}

// Equations refer to Pinski et al. (JCP 143, 034108, 2015; DOI: 10.1063/1.4926879)

class DLPNO : public Wavefunction {
   protected:
    /// what quantum chemistry module are we running
    DLPNOMethod algorithm_;
    
    /// threshold for PAO domain size
    double T_CUT_DO_;
    /// threshold for PNO truncation
    double T_CUT_PNO_;
    /// trace threshold for PNO truncation (CC only)
    double T_CUT_TRACE_;
    /// pair energy threshold for PNO truncation (CC only)
    double T_CUT_ENERGY_;
    /// threshold for PNO truncation for MP2 pairs (for DLPNO-CC methods)
    double T_CUT_PNO_MP2_;
    /// trace threshold for PNO truncation for MP2 pairs (for DLPNO-CC methods)
    double T_CUT_TRACE_MP2_;
    /// pair energy threshold for PNO truncation for MP2 pairs (for DLPNO-CC methods)
    double T_CUT_ENERGY_MP2_;
    /// tolerance to separate pairs into CCSD and MP2 pairs
    double T_CUT_PAIRS_;
    /// tolerance to separate MP2 pairs in between crude and refined prescreening
    double T_CUT_PAIRS_MP2_;
    /// tolerance for energy of a pair for it to be considered a "dipole pair"
    double T_CUT_PRE_;
    /// tolerance for local density fitting (by Mulliken population)
    double T_CUT_MKN_;
    /// T_CUT_PNO scaling factor for diagonal PNOs
    double T_CUT_PNO_DIAG_SCALE_;
    /// T_CUT_PNO scaling for core orbitals
    double T_CUT_PNO_CORE_SCALE_;
    /// Tolerance for TNO truncation for triples (by occupation number)
    double T_CUT_TNO_;
    
    /// toggle core and disk options based on available memory
    bool toggle_memory_;
    /// number of core orbitals (0 if freeze_core = True)
    int ncore_;

    /// auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;
    SharedMatrix full_metric_;
    std::vector<double> J_metric_shell_diag_; ///< used in AO ERI screening

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
    SharedMatrix DOI_uv_; // PAO/PAO

    // approximate LMO/LMO pair energies from dipole integrals (EQ 17)
    // used to screen out and estimate weakly interacting LMO/LMO pairs
    SharedMatrix dipole_pair_e_; ///< actual approximate pair energy (used in final energy calculation)
    SharedMatrix dipole_pair_e_bound_; ///< upper bound to approximate pair energies (used for screening)

    /// How much memory is used by storing each of the DF integral types
    size_t qij_memory_;
    size_t qia_memory_;
    size_t qab_memory_;

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
    std::vector<double> occ_pno_;       ///< lowest PNO occupation number per PNO
    std::vector<double> trace_pno_;     ///< total trace(Dij) recovered per PNO
    std::vector<double> e_ratio_pno_;   ///< percentage of correlation energy recovered by PNOs
    std::vector<double> de_pno_;   ///< PNO truncation energy error
    std::vector<double> de_pno_os_;   ///< opposite-spin contributions to de_pno_
    std::vector<double> de_pno_ss_;   ///< same-spin contributions to de_pno_

    /// pre-screening energies
    double de_dipole_; ///< energy correction for distant (LMO, LMO) pairs
    double de_pno_total_; ///< energy correction for PNO truncation (total)
    double de_pno_total_os_; ///< energy correction for PNO truncation (opposite-spin)
    double de_pno_total_ss_; ///< energy correction for PNO truncation (same-spin)
    double e_lmp2_non_trunc_; ///< LMP2 energy in a pure PAO basis (Strong and Weak Pairs Only)
    double e_lmp2_trunc_; ///< LMP2 energy computed with (truncated) PNOs (Strong Pairs Only)
    double de_lmp2_eliminated_; ///< LMP2 correction for eliminated pairs (surviving pairs after dipole screening that
    // are neither weak nor strong)
    double de_weak_; ///< Energy contribution for weak pairs

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

    /* Takes an atom index and a global LMO index 
    to return the sparse LMO index on that atom,
    (-1) if that LMO is not on the riatom's extended domain */
    std::vector<std::vector<int>> riatom_to_lmos_ext_dense_;
    /* Takes an atom index and a global PAO index 
    to return the sparse PAO index on that atom,
    (-1) if that PAO is not on the riatom's extended domain */
    std::vector<std::vector<int>> riatom_to_paos_ext_dense_;
    std::vector<std::vector<bool>> riatom_to_atoms1_dense_;
    std::vector<std::vector<bool>> riatom_to_atoms2_dense_;
    std::vector<std::vector<int>> lmopair_to_lmos_dense_;

    /// Useful for generating DF integrals (TODO: Replace this with "index_list" function)
    std::vector<std::vector<std::vector<int>>> lmopair_lmo_to_riatom_lmo_;
    std::vector<std::vector<std::vector<int>>> lmopair_pao_to_riatom_pao_;
    std::vector<std::vector<std::pair<int,int>>> riatom_to_pao_pairs_; ///< Which (u,v) pao pairs belong to an riatom
    std::vector<std::vector<std::vector<int>>> riatom_to_pao_pairs_dense_; ///< For each riatom, returns the index of the element in qab tensor

    /// PSIO object (helps with reading/writing large tensors)
    std::shared_ptr<PSIO> psio_;

    void common_init();

    // Helper functions
    void C_DGESV_wrapper(SharedMatrix A, SharedMatrix B);

    std::pair<SharedMatrix, SharedVector> canonicalizer(SharedMatrix C, SharedMatrix F);
    std::pair<SharedMatrix, SharedVector> orthocanonicalizer(SharedMatrix S, SharedMatrix F);

    SharedVector flatten_mats(const std::vector<SharedMatrix>& mat_list);

    void copy_flat_mats(SharedVector flat, std::vector<SharedMatrix>& mat_list);

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
    void prep_sparsity(bool initial, bool final);

    /// Compute the auxiliary metric
    void compute_metric();
    /// Compute three-index integrals in LMO/LMO basis with linear scaling
    void compute_qij();
    /// Compute three-index integrals in LMO/PAO basis with linear scaling (EQ 11)
    void compute_qia();
    /// Compute three-index integrals in PAO/PAO basis with linear scaling
    void compute_qab();

    /// form pair exch operators (EQ 15) and SC amplitudes (EQ 18); transform to PNO basis
    void pno_transform();

    void print_aux_domains();
    void print_pao_domains();
    void print_lmo_domains(bool initial);
    void print_aux_pair_domains();
    void print_lmo_pair_domains();
    void print_pao_pair_domains();

   public:
    DLPNO(SharedWavefunction ref_wfn, Options& options);
    ~DLPNO() override;

    double compute_energy() override;
};

// Equations refer to Pinski et al. (JCP 143, 034108, 2015; DOI: 10.1063/1.4926879)

class DLPNOMP2 : public DLPNO {
   protected:
    // PNO overlap matrices
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik_; ///< pnooverlaps
    
    // final energies
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

// Equations refer to Jiang et al. (JCP 161, 082502, 2024; DOI: 10.1063/5.0219963)
class PSI_API DLPNOCCSD : public DLPNO {
   protected:
    /// Use low memory algorithm to store PNO overlaps?
    bool low_memory_overlap_;
    /// Write (Q_ij | m_ij a_ij) integrals to disk?
    bool write_qia_pno_;
    /// Write (Q_ij | a_ij b_ij) integrals to disk?
    bool write_qab_pno_;

    /// PNO overlap integrals
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_nn_; ///< pno overlaps
    std::vector<std::vector<SharedMatrix>> S_pno_ij_mn_; ///< pno overlaps

    /// Coupled-cluster amplitudes
    std::vector<SharedMatrix> T_ia_; ///< singles amplitudes [naocc x (npno_ii, 1)]
    std::vector<SharedMatrix> T_n_ij_; ///< projected singles amplitudes [n_lmo_pairs x (nlmo_ij, npno_ij)] (Jiang Eq. 70)

    // => Strong and Weak Pair Info <=//

    std::vector<double> e_ij_mp2_scale_; ///< how much to scale MP2 energies for scaled approximation to PAO-LMP2

    std::vector<std::vector<int>> i_j_to_ij_strong_; ///< LMO indices (i, j) to significant strong pair index (ij); insignificant (i, j) maps to -1
    std::vector<std::pair<int,int>> ij_to_i_j_strong_;  ///< LMO strong pair index (ij) to both LMO indices (i, j)
    std::vector<int> ij_to_ji_strong_; ///< LMO strong pair index (ij) to LMO pair index (ji)

    std::vector<std::vector<int>> i_j_to_ij_weak_; ///< LMO indices (i, j) to significant weak pair index (ij); insignificant (i, j) maps to -1
    std::vector<std::pair<int,int>> ij_to_i_j_weak_;  ///< LMO weak pair index (ij) to both LMO indices (i, j)
    std::vector<int> ij_to_ji_weak_; ///< LMO weak pair index (ij) to LMO pair index (ji)

    // => CCSD Integrals <= //

    // Uses the notation of Jiang et al. Eq. 71-74
    // J_pqrs_ => J_{pq}^{rs} = (pq|rs), K_prqs_ => K_{pq}^{rs} = (pr|qs)
    // L_prqs_ => L_{pq}^{rs} = 2(pr|qs) - (ps|qr), M_{pq}^{rs} = 2(pr|qs) - (pq|rs)
    
    /// 1-external integrals
    std::vector<SharedMatrix> K_mibj_; /// (m_{ij} i | b_{ij} j)
    std::vector<SharedMatrix> J_ijmb_; /// (i j | m_{ij} b_{ij})
    std::vector<SharedMatrix> L_mibj_; /// 2.0 (m_{ij} i | b_{ij} j) - (m_{ij} j | b_{ij} i)

    /// 2-external integrals
    std::vector<SharedMatrix> L_iajb_; /// 2.0 * (i a_{ij} | j b_{ij}) - (i b_{ij} | j a_{ij})

    /// 2-external non-projected integrals
    std::vector<std::vector<SharedMatrix>> J_ikac_non_proj_; /// (i k | a_{ij} c_{kj})
    std::vector<std::vector<SharedMatrix>> K_iakc_non_proj_; /// (i a_{ij} | k c_{kj})

    /// 3-external integrals
    std::vector<SharedMatrix> K_ivvv_; /// (i e_{ij} | a_{ij} f_{ij}) (stored as (e, a * f))

    // Density-fitted integrals (only computed over strong pairs)
    std::vector<std::vector<SharedMatrix>> Qma_ij_; // (Q_{ij} | m_{ij} a_{ij})
    std::vector<std::vector<SharedMatrix>> Qab_ij_; // (Q_{ij} | a_{ij} b_{ij})

    std::vector<SharedMatrix> i_Qk_ij_;   // (Q_{ij} | k_{ij} i)
    std::vector<SharedMatrix> i_Qa_ij_;   // (Q_{ij} | a_{ij} i)
    std::vector<SharedMatrix> i_Qk_t1_;   // (Q_{ij} | k_{ij} i) [T1-dressed] (Jiang Eq. 91)
    std::vector<SharedMatrix> i_Qa_t1_;   // (Q_{ij} | a_{ij} i) [T1-dressed] (Jiang Eq. 92)

    // Dressed Fock matrices (used in DLPNO-T1-CCSD)
    SharedMatrix Fkj_; // Jiang Eq. 94
    std::vector<SharedMatrix> Fkc_; // Jiang Eq. 95
    std::vector<SharedMatrix> Fai_; // Jiang Eq. 96
    std::vector<SharedMatrix> Fab_; // Jiang Eq. 97

    double e_lmp2_; ///< raw (uncorrected) local MP2 correlation energy
    double e_lccsd_; ///< raw (uncorrected) local CCSD correlation energy
    /// Baseline and temporary-memory estimates retained from DLPNOCCSD::estimate_memory(), in doubles
    size_t ccsd_baseline_memory_doubles_ = 0;
    size_t ccsd_iteration_workspace_doubles_ = 0;
    size_t ccsd_peak_memory_doubles_ = 0;

    /// Returns the appropriate overlap matrix given two LMO pairs
    inline SharedMatrix S_PNO(const int ij, const int mn);
    /// Encapsulates the reading in of (Q_{ij}|m_{ij} a_{ij}) integrals (regardless of core or disk)
    inline std::vector<SharedMatrix> QIA_PNO(const int ij);
    /// Encapsulates the reading in of (Q_{ij}|a_{ij} b_{ij}) integrals (regardless of core or disk)
    inline std::vector<SharedMatrix> QAB_PNO(const int ij);

    /// These functions split up pairs that survive the initial dipole screening
    // The "crude" pre-screening step splits up semi-canonical MP2 pairs from the rest,
    // while the non-crude splits up strong/weak pairs (see Riplinger 2016 for more details)
    template<bool crude> void pair_prescreening();
    template<bool crude> std::vector<double> compute_pair_energies();
    template<bool crude> double filter_pairs(const std::vector<double>& e_ijs);

    /// Runs preceeding DLPNO-MP2 computation before DLPNO-CCSD iterations
    std::vector<double> pno_lmp2_iterations();
    /// Recompute PNOs after DLPNO-MP2 converges
    void recompute_pnos();

    /// compute PNO/PNO overlap matrices for DLPNO-CCSD
    void compute_pno_overlaps();

    // => Computing integrals <= //

    /// A function to estimate integral memory costs
    void estimate_memory();
    /// Compute PNO integrals for CC computations
    void compute_pno_integrals();

    // => CCSD intermediates <= //

    /// Jiang Equation 82
    std::vector<SharedMatrix> compute_beta();
    /// Jiang Equation 83
    std::vector<SharedMatrix> compute_gamma();
    /// Jiang Equation 84
    std::vector<SharedMatrix> compute_delta();
    /// Jiang Equation 86
    SharedMatrix compute_Fkj_double_tilde();

    /// compute T1-dressed DF integrals (Jiang Eq. 91-92)
    void t1_ints();
    /// compute T1-dressed Fock matrix intermediates (Jiang Eq. 94-101)
    void t1_fock();

    /// computes singles residuals in LCCSD equations, using pre-allocated memory (Jiang Eq. 32)
    void compute_R_ia(std::vector<SharedMatrix>& R_ia, std::vector<std::vector<SharedMatrix>>& R_ia_buffer);
    /// computes doubles residuals in LCCSD equations, using pre-allocated memory (Jiang Eq. 19)
    void compute_R_iajb(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb);

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

// Equations refer to Jiang et al. (JCP 161, 082502, 2024; DOI: 10.1063/5.0219963)

class PSI_API DLPNOCCSD_T : public DLPNOCCSD {
   protected:
    // Sparsity information, NOTE: only unique triplets i <= j <= k are used
    SparseMap lmotriplet_to_ribfs_; ///< which ribfs are on an LMO triplet (i, j, k)
    SparseMap lmotriplet_to_lmos_; ///< which LMOs l form a significant pair with (i, j, or k)
    SparseMap lmotriplet_to_paos_; ///< which PAOs span the virtual space of a triplet of LMOs?
    std::unordered_map<size_t, int> i_j_k_to_ijk_; ///< LMO indices (i, j, k) to significant LMO triplet index (ijk), -1 if not found
    std::vector<std::tuple<int, int, int>> ijk_to_i_j_k_; ///< LMO triplet index (ijk) to LMO index tuple (i, j, k)

    /// triplet natural orbitals (TNOs)
    std::vector<SharedMatrix> W_iajbkc_; ///< W3 intermediate for each lmo triplet (Jiang Eq. 109)
    std::vector<SharedMatrix> V_iajbkc_; ///< V3 intermeidate for each lmo triplet (Jiang Eq. 110)
    std::vector<SharedMatrix> T_iajbkc_; ///< Triples amplitude for each lmo triplet
    std::vector<SharedMatrix> X_tno_; ///< global PAO -> canonical TNO transforms
    std::vector<SharedVector> e_tno_; ///< TNO orbital energies
    std::vector<int> n_tno_; ///< number of tnos per triplet domain
    std::vector<double> e_ijk_; ///< energy of triplet ijk (used for pre-screening and convergence purposes)
    std::vector<double> tno_cutoff_; ///< absolute TNO cutoff for each energy-ranked triplet in iterative (T)
    std::vector<bool> is_strong_triplet_; ///< whether or not triplet is strong

    /// Write intermediates (W and V) to disk?
    bool write_intermediates_ = false;
    /// Write amplitudes to disk?
    bool write_amplitudes_ = false;

    /// final energies
    double de_lccsd_t_screened_; ///< energy contribution from screened triplets
    double e_lccsd_t_; ///< local (T) correlation energy
    double E_T_; ///< raw iterative (T) energy at weaker triples cutoffs

    /// Create sparsity maps for triples
    void triples_sparsity(bool prescreening);
    /// Create TNOs (Triplet Natural Orbitals) for DLPNO-(T)
    void tno_transform(double tno_tolerance, bool use_tuple_cutoffs = false);
    /// Sort triplets to split between "strong" and "weak" triplets (for (T) iterations)
    void sort_triplets(double e_total);

    /// A helper function to transform triples-like tensor
    SharedMatrix matmul_3d(SharedMatrix A, SharedMatrix X, int dim_old, int dim_new);
    /// Returns a symmetrized version of a triples-like tensor (in i <= j <= k ordering)
    SharedMatrix triples_permuter(const SharedMatrix& X, int i, int j, int k, bool reverse=false);
    /// compute (T) iteration energy (Jiang Eq. 53)
    double compute_t_iteration_energy();

    /// L_CCSD(T0) energy (Jiang Eq. 53, 109-110)
    double compute_lccsd_t0(bool save_memory=false);
    /// A function to estimate (T) memory costs
    void estimate_memory();
    /// L_CCSD(T) iterations (Jiang Eq. 111-112)
    double lccsd_t_iterations();

    void print_header();

    void print_results();

   public:
    DLPNOCCSD_T(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSD_T() override;

    double compute_energy() override;
};

#ifdef USING_Einsums

// Equations for DLPNOCCSDT refer to Jiang et al. (JCTC 21, 2386, 2025;
// DOI: 10.1021/acs.jctc.4c01716) and its Supporting Information unless noted otherwise.
class DLPNOCCSDT : public DLPNOCCSD_T {
   protected:
    // DF integrals in the domain of triplet ijk
    std::vector<einsums::Tensor<double, 2>> q_io_;
    std::vector<einsums::Tensor<double, 2>> q_jo_;
    std::vector<einsums::Tensor<double, 2>> q_ko_;

    std::vector<einsums::Tensor<double, 2>> q_iv_;
    std::vector<einsums::Tensor<double, 2>> q_jv_;
    std::vector<einsums::Tensor<double, 2>> q_kv_;

    std::vector<einsums::Tensor<double, 3>> q_ov_;
    std::vector<einsums::Tensor<double, 3>> q_vv_;

    // Write expensive integrals (q_{ijk}| m_{ijk} a_{ijk}) and (q_{ijk} | a_{ijk} b_{ijk}) to disk?
    bool disk_ints_;
    // Damping ratio (how much of the original triples amplitude to keep)
    double damping_ratio_;

    // Singles Amplitudes
    std::vector<einsums::Tensor<double, 2>> T_n_ijk_;
    // Einsums clone of Psi4 T3 amplitudes
    std::vector<einsums::Tensor<double, 3>> T_iajbkc_clone_;
    // Contravariant (spin-summed) triples amplitudes; manuscript Eq. (79)
    std::vector<einsums::Tensor<double, 3>> U_iajbkc_;

    /// Cached (ia|jb), (ia|kb), and (ja|kb) TNO-domain exchange matrices.
    /// They are invariant during the CCSDT iterations and are shared by the
    /// separate T3 -> R1 and T3 -> R2 contraction passes.
    std::vector<std::array<einsums::Tensor<double, 2>, 3>> K_ovov_tno_cache_;

    // List of triples sorted by number of TNOs
    std::vector<int> sorted_triplets_;
    // Number of threads
    int nthread_;
    // Energy expression
    double e_lccsdt_;

    /// Persistent in-core state at the end of the CCSDT stage, in doubles.
    /// The CCSDT(Q) estimator uses this when CCSDTQ must retain all lower-rank data.
    size_t ccsdt_resident_memory_doubles_ = 0;
    /// Largest temporary lower-rank residual workspace required by a CCSDT iteration, in doubles.
    /// The CCSDTQ estimator combines this with its T4-dependent residual workspaces.
    size_t ccsdt_iteration_workspace_doubles_ = 0;

    /// Load disk-backed (Q_{ijk}|m_{ijk} a_{ijk}) integrals into q_ov_[ijk].
    void load_qia_tno(int ijk);
    /// Load disk-backed (Q_{ijk}|a_{ijk} b_{ijk}) integrals into q_vv_[ijk].
    void load_qab_tno(int ijk);

    /// Helper function for transforming amplitudes from one TNO space to another
    einsums::Tensor<double, 3> matmul_3d_einsums(const einsums::Tensor<double, 3> &A, const SharedMatrix &X, int dim_old, int dim_new);
    /// Helper function for transforming amplitudes from one TNO space to another (one index only)
    einsums::Tensor<double, 3> matmul_3d_index(const einsums::Tensor<double, 3> &A, const SharedMatrix &X, int index);
    /// Apply the direct or conjugate triples permutation operators, manuscript Eqs. (81)-(86) or (91)-(96)
    einsums::Tensor<double, 3> triples_permuter_einsums(const einsums::Tensor<double, 3> &X, int i, int j, int k, bool reverse=false);
    /// Spin-sum a triples orbital tensor following Matthews and Stanton, JCP 142, 064108 (2015)
    einsums::Tensor<double, 3> triples_spin_summation(const einsums::Tensor<double, 3> &X);
    /// Recover an equivalent orbital triples tensor using Matthews and Stanton Eq. (27)
    einsums::Tensor<double, 3> triples_spin_desummation(const einsums::Tensor<double, 3> &X);

    /// Factorized cross-domain contribution implied by rho_k(d,b,c), manuscript Eq. (101) and SI Algorithm S3
    std::vector<einsums::Tensor<double, 3>> rho_dbck_contribution();
    /// Add triples to the singles residual: canonical Eq. (32), DLPNO Eq. (80), Algorithm 1
    void compute_R_ia_triples(std::vector<SharedMatrix>& R_ia, std::vector<std::vector<SharedMatrix>>& R_ia_buffer);
    /// Add triples to the doubles residual: canonical Eqs. (33)-(34), DLPNO Eqs. (87)-(90), Algorithm 2
    void compute_R_iajb_triples(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb, std::vector<std::vector<SharedMatrix>>& R_iajb_buffer);
    /// Form the triples residual: canonical Eqs. (37)-(50), DLPNO Eqs. (99)-(110), Algorithms 3-4
    void compute_R_iajbkc(std::vector<SharedMatrix>& R_iajbkc);

    void print_header();
    void estimate_memory();
    void compute_integrals();
    void lccsdt_iterations();
    void print_results();

   public:
    DLPNOCCSDT(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSDT() override;

    double compute_energy() override;
};

// Equations and algorithms for DLPNOCCSDT_Q refer to Jiang, Schaefer, and Turney,
// J. Chem. Phys. 162, 144102 (2025), DOI: 10.1063/5.0257672.
class DLPNOCCSDT_Q : public DLPNOCCSDT {
   protected:
    struct QuadrupletEnergyIntermediates {
        /// g_i'(a,b,e) = (i'a|be), Algorithm 1; used by Eq. (19), term 1, and Algorithm 4.
        std::array<einsums::Tensor<double, 3>, 4> K_iabe;
        /// g_i'j'(a,m) = (i'a|j'm), Algorithm 1; used by Eq. (19), term 2, and Algorithm 4.
        std::array<einsums::Tensor<double, 2>, 16> K_iajm;
        /// g_i'j'(a,b) = (i'a|j'b), Algorithm 1; used by canonical Eqs. (25)-(26).
        std::array<einsums::Tensor<double, 2>, 10> K_iajb;
        /// Projected contravariant doubles U_i'j'(a,b), Algorithm 1 and canonical Eq. (25).
        std::array<einsums::Tensor<double, 2>, 10> U_iajb;
    };

    // All 24 permutations of four occupied-orbital positions.
    constexpr static std::array<std::tuple<int, int, int, int>, 24> quadruple_permutations_ = {std::make_tuple(0, 1, 2, 3), std::make_tuple(0, 1, 3, 2), 
        std::make_tuple(0, 2, 1, 3), std::make_tuple(0, 2, 3, 1), std::make_tuple(0, 3, 1, 2), std::make_tuple(0, 3, 2, 1), 
        std::make_tuple(1, 0, 2, 3), std::make_tuple(1, 0, 3, 2), std::make_tuple(1, 2, 0, 3), std::make_tuple(1, 2, 3, 0), 
        std::make_tuple(1, 3, 0, 2), std::make_tuple(1, 3, 2, 0), std::make_tuple(2, 0, 1, 3), std::make_tuple(2, 0, 3, 1), 
        std::make_tuple(2, 1, 0, 3), std::make_tuple(2, 1, 3, 0), std::make_tuple(2, 3, 0, 1), std::make_tuple(2, 3, 1, 0), 
        std::make_tuple(3, 0, 1, 2), std::make_tuple(3, 0, 2, 1), std::make_tuple(3, 1, 0, 2), std::make_tuple(3, 1, 2, 0), 
        std::make_tuple(3, 2, 0, 1), std::make_tuple(3, 2, 1, 0)};

    // Packed upper triangle of the four occupied positions. Pair tensors use
    // X_ij(a,b) = X_ji(b,a), so ten components replace the former 4 x 4 arrays.
    constexpr static std::array<std::array<int, 4>, 4> pair_positions_ = {
        std::array<int, 4>{0, 1, 2, 3}, std::array<int, 4>{1, 4, 5, 6},
        std::array<int, 4>{2, 5, 7, 8}, std::array<int, 4>{3, 6, 8, 9}};
    constexpr static int pair_position(int i, int j) { return pair_positions_[i][j]; }

    SparseMap lmoquadruplet_to_ribfs_; ///< RI basis functions in each LMO quadruplet domain (i, j, k, l)
    SparseMap lmoquadruplet_to_lmos_; ///< which LMOs m form a significant pair with (i, j, k, or l)
    SparseMap lmoquadruplet_to_paos_; ///< which PAOs span the virtual space of a quadruplet of LMOs?
    std::unordered_map<size_t, int> i_j_k_l_to_ijkl_; ///< LMO indices (i, j, k, l) to significant LMO quadruplet index (ijkl), -1 if not found
    std::vector<std::tuple<int, int, int, int>> ijkl_to_i_j_k_l_; ///< LMO quadruplet index (ijkl) to LMO index tuple (i, j, k, l)
    std::vector<std::tuple<int, int, int, int>> ijkl_to_i_j_k_l_full_; ///< LMO quadruplet indices with no i <= j <= k <= l restriction
    std::vector<int> sorted_quadruplets_; ///< quadruplets sorted by number of QNOs

    /// Quadruples natural orbitals (QNOs), formed from the six-pair density of Eq. (41).
    std::vector<einsums::Tensor<double, 4>> T_iajbkcld_; ///< T4 amplitudes, canonical Eqs. (20), (23)-(24)
    std::vector<einsums::Tensor<double, 4>> gamma_ijkl_; ///< T4 source term, canonical Eq. (19), DLPNO Algorithm 2
    /// Reusable inputs to the [Q] and (Q) energy contractions, grouped by quadruplet
    /// so their lifetime and optional disk backing can be managed as one unit.
    std::vector<QuadrupletEnergyIntermediates> quadruplet_energy_intermediates_;
    std::vector<SharedMatrix> X_qno_; ///< PAO -> canonical QNO transforms
    std::vector<SharedVector> e_qno_; ///< QNO orbital energies
    std::vector<int> n_qno_; ///< number of qnos per quadruplet domain
    std::vector<double> e_ijkl_; ///< energy of quadruplet ijkl (used for pre-screening and convergence purposes)
    std::vector<double> ijkl_scale_; ///< Scaling factor applied to quadruplet energy ijkl based on MP2 scaling
    std::vector<double> qno_cutoff_; ///< absolute QNO cutoff for each energy-ranked quadruplet in iterative (Q)
    std::vector<bool> is_strong_quadruplet_; ///< whether or not quadruplet is strong

    /// Store iterative-(Q) source, amplitudes, and reusable energy intermediates on disk.
    bool write_quadruples_intermediates_ = false;

    /// final energies
    double de_lccsdt_q_screened_; ///< energy contribution from screened quadruples
    double e_lccsdt_q_; ///< local (Q) correlation energy
    double E_Q_; ///< raw iterative (Q) energy at weaker quadruples cutoffs

    /// Create sparsity maps for quadruples
    void quadruples_sparsity(bool prescreening);
    /// Create QNOs (Quadruplet Natural Orbitals) for DLPNO-(Q)
    void qno_transform(double qno_tolerance, bool use_tuple_cutoffs = false);
    /// Sort quadruplets to split between "strong" and "weak" quadruplets (for (Q) iterations)
    void sort_quadruplets();

    /// Transform all four virtual indices between QNO spaces (semidirect overlaps, Eqs. (58)-(61)).
    einsums::Tensor<double, 4> matmul_4d(const einsums::Tensor<double, 4>& A,
                                         const SharedMatrix& X, int dim_old, int dim_new);
    /// Transform in canonical storage order and apply the occupied-column
    /// permutation only in the smaller target space, avoiding a permuted
    /// dim_old^4 input copy while keeping all four contractions as GEMMs.
    einsums::Tensor<double, 4> matmul_4d_permuted(
        const einsums::Tensor<double, 4>& A, const SharedMatrix& X, int dim_old,
        int dim_new, int i, int j, int k, int l);
    /// Apply the occupied/QNO permutation maps of Eqs. (51)-(55).
    einsums::Tensor<double, 4> quadruples_permuter(const einsums::Tensor<double, 4>& X, int i, int j, int k, int l);

    /// Save/load rank-four tensors and the reusable Algorithm 1 energy bundle through PSIO.
    void save_quadruples_tensor(const einsums::Tensor<double, 4>& tensor, const std::string& label, int ijkl);
    einsums::Tensor<double, 4> load_quadruples_tensor(const std::string& label, int ijkl);
    void save_quadruplet_energy_intermediates(const QuadrupletEnergyIntermediates& intermediates, int ijkl);
    QuadrupletEnergyIntermediates load_quadruplet_energy_intermediates(int ijkl);

    /// Form Gamma and semicanonical T4 (canonical Eqs. (19)-(20), DLPNO Algorithms 1-2)
    /// and return the semicanonical (Q0) energy of Eqs. (25)-(26).
    double compute_gamma_ijkl(bool store_amplitudes=false);
    /// Evaluate the [Q] and (Q) energy contractions (canonical Eqs. (25)-(26), Algorithms 3-4).
    double compute_quadruplet_energy(int ijkl, const einsums::Tensor<double, 4>& T4,
                                     const QuadrupletEnergyIntermediates& intermediates);
    /// Estimate iterative-(Q) resident and thread-workspace peaks and select disk backing.
    void estimate_memory();
    /// Iterate the non-semicanonical T4 equations (canonical Eqs. (23)-(24), DLPNO Algorithm 5).
    double lccsdt_q_iterations();

    void print_header();
    void print_results();

   public:
    DLPNOCCSDT_Q(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSDT_Q() override;

    double compute_energy() override;
};

// Equations for DLPNOCCSDTQ refer to Jiang et al. (JCTC 22, 2825--2845, 2026;
// DOI: 10.1021/acs.jctc.5c01910) and its Supporting Information unless noted otherwise.
class DLPNOCCSDTQ : public DLPNOCCSDT_Q {
   protected:
    // DF domain integrals
    std::vector<std::array<einsums::Tensor<double, 2>, 4>> q_io_list_; ///< (Q_{ijkl} | [i, j, k, l] m_{ijkl})
    std::vector<std::array<einsums::Tensor<double, 2>, 4>> q_iv_list_; ///< (Q_{ijkl} | [i, j, k, l] a_{ijkl})
    std::vector<einsums::Tensor<double, 3>> q_ov_ijkl_; ///< (Q_{ijkl} | m_{ijkl} a_{ijkl})
    std::vector<einsums::Tensor<double, 3>> q_vv_ijkl_; ///< (Q_{ijkl} | a_{ijkl} b_{ijkl})

    // Extended PNO (XPNO) information
    SparseMap lmopair_to_paos_ext_;           ///< LMO pair to extended PAO domain
    std::vector<SharedMatrix> X_xpno_;         ///< Extended-PAO-to-canonical-XPNO transforms
    std::vector<SharedVector> e_xpno_;         ///< Canonical XPNO orbital energies
    std::vector<int> n_xpno_;                  ///< Number of XPNOs in each extended pair domain

    /// Load disk-backed (Q_{ijkl} | m_{ijkl} a_{ijkl}) integrals into q_ov_ijkl_[ijkl].
    void load_qia_qno(int ijkl);
    /// Load disk-backed (Q_{ijkl} | a_{ijkl} b_{ijkl}) integrals into q_vv_ijkl_[ijkl].
    void load_qab_qno(int ijkl);

    /// Write the largest QNO-basis integral tensors to disk (enabled by default).
    bool disk_qno_integrals_ = true;
    /// Include T4 amplitudes and residuals in the on-disk DIIS vectors.
    bool extrapolate_t4_ = true;
    /// Recompute one XPNO-projected T4 block at a time instead of retaining the
    /// complete [kl][mn] bank when the memory estimator selects low-memory mode.
    bool stream_xpno_t4_ = false;
    /// Fraction of the preceding T4 amplitude retained when damping is activated.
    double quadruples_damping_ratio_ = 0.0;
    /// Local CCSDTQ correlation energy.
    double e_lccsdtq_ = 0.0;

    /// Singles amplitudes t_m^a projected into the QNO space of ijkl.
    std::vector<einsums::Tensor<double, 2>> T_m_ijkl_;
    /// T_mnkl projected into the XPNO space of kl (main-text Eq. (83)).
    std::vector<std::vector<einsums::Tensor<double, 4>>> T_mnkl_xpno_;

    /// Form alpha_ijkl from T_ijkl (main-text Eq. (30)).
    inline einsums::Tensor<double, 4> form_alpha_ijkl(const einsums::Tensor<double, 4>& T_ijkl);
    /// Form beta_ijkl from alpha_ijkl (main-text Eq. (31)).
    inline einsums::Tensor<double, 4> form_beta_ijkl(const einsums::Tensor<double, 4>& alpha_ijkl);
    /// Form the spin-summed quadruples tensor used by the closed-shell equations.
    einsums::Tensor<double, 4> quadruples_spin_summation(const einsums::Tensor<double, 4> &X);
    /// Apply the minimum-norm spin pseudoinverse of Matthews and Stanton, Eq. (28)
    /// (JCP 142, 064108, 2015; DOI: 10.1063/1.4907278), removing the redundant
    /// components of the nonorthogonal spin-adapted quadruples representation.
    einsums::Tensor<double, 4> quadruples_spin_desummation(const einsums::Tensor<double, 4> &X);

    /// Flatten Psi4 matrices and, optionally, native Einsums T4/R4 tensors into one DIIS vector.
    SharedVector flatten_ccsdtq_diis(const std::vector<SharedMatrix>& matrices,
                                     const std::vector<einsums::Tensor<double, 4>>& rank4_tensors,
                                     bool include_t4) const;
    /// Scatter an extrapolated DIIS vector back into the Psi4 matrices and native T4/R4 tensors.
    void copy_ccsdtq_diis(const SharedVector& flat, std::vector<SharedMatrix>& matrices,
                          std::vector<einsums::Tensor<double, 4>>& rank4_tensors, bool include_t4) const;

    /// Create extended pair natural orbitals (XPNOs) for the T_mnkl contractions.
    void xpno_transform(double xpno_tolerance);

    /// Form projected T_{mnkl}^{abcd} amplitudes in the XPNO domain of kl
    /// to reduce the cost of the computationally dominant contractions.
    void form_T_mnkl_xpno();

    /// Form the pair-domain delta-L and delta-M corrections together, sharing
    /// streamed (me|nf) slices (main-text Eqs. (98)--(102); SI Eqs. (S13),(S15)).
    std::pair<std::vector<einsums::Tensor<double, 4>>,
              std::vector<einsums::Tensor<double, 4>>>
    compute_delta_L_M();

    /// Add the T4-dependent doubles residual (main-text Eq. (93); SI Eq. (S16)).
    void add_t4_to_doubles_residual(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb,
                                    std::vector<std::vector<SharedMatrix>>& R_iajb_buffer);
    /// Add the T4-dependent triples residual (main-text Eqs. (94)--(96); SI Eqs. (S17)--(S19)).
    void add_t4_to_triples_residual(std::vector<SharedMatrix>& R_iajbkc);
    /// Form the local quadruples residual (main-text Eq. (50); SI Eqs. (S20)--(S34)).
    void compute_quadruples_residual(std::vector<einsums::Tensor<double, 4>>& R_iajbkcld);

    /// Print full-quadruples thresholds, iteration controls, and requested storage policy.
    void print_header();
    /// Estimate resident and per-thread peaks and select safe automatic storage fallbacks.
    void estimate_memory();
    /// Build the QNO-basis density-fitted integrals required by the CCSDTQ iteration.
    void compute_integrals();
    /// Iterate the coupled local T1/T2/T3/T4 amplitude equations.
    void lccsdtq_iterations();
    /// Print the final energy decomposition and post-CCSDT increments.
    void print_results();

   public:
    DLPNOCCSDTQ(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOCCSDTQ() override;

    double compute_energy() override;
};

#endif  // USING_Einsums

}
}

#endif
