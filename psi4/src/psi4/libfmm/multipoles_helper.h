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

#ifndef libfmm_fmm_tree_H
#define libfmm_fmm_tree_H

#include "psi4/pragma.h"

#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libfmm/multipoles_helper.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <unordered_map>

// Extent normalization constant (to calculate extents of shell pairs)
#define ENC (1.306076308436251)

namespace psi {

class Options;

enum class ContractionType {
    DIRECT,
    DF_AUX_PRI,
    DF_PRI_AUX,
    METRIC
};

class PSI_API ShellPair {
    protected:
      // The basis set of the first shell-pair index
      std::shared_ptr<BasisSet> bs1_;
      // The basis set of the second shell-pair index
      std::shared_ptr<BasisSet> bs2_;
      // The index of the shell-pair
      std::pair<int, int> pair_index_;
      // Exponent of most diffuse basis function in shell pair
      double exp_;
      // Center of shell pair (As defined in bagel FMM as the average)
      Vector3 center_;
      // Radial extent of shell pair
      double extent_;
      // The multipole moments (per basis pair (pq) the shell pair (PQ)), centered at the lowest level box the shell belongs to
      std::vector<std::shared_ptr<RealSolidHarmonics>> mpoles_;
      // Multipole coefficients of shellpair
      std::shared_ptr<HarmonicCoefficients> mpole_coefs_;

    public:
      ShellPair(std::shared_ptr<BasisSet>& bs1, std::shared_ptr<BasisSet>& bs2, std::pair<int, int> pair_index, 
                std::shared_ptr<HarmonicCoefficients>& mpole_coefs, double cfmm_extent_tol);

      // Calculate the multipole moments of the Shell-Pair about a center
      void calculate_mpoles(Vector3 box_center, std::shared_ptr<OneBodyAOInt>& s_ints,
                            std::shared_ptr<OneBodyAOInt>& mpole_ints, int lmax);

      // Returns the shell pair index
      std::pair<int, int> get_shell_pair_index() { return pair_index_; }
      // Returns the center of the shell pair
      Vector3 get_center() { return center_; }
      // Returns the radial extent of the shell pair
      double get_extent() { return extent_; }
      // Returns the multipole moments of the shell pairs about a center
      std::vector<std::shared_ptr<RealSolidHarmonics>>& get_mpoles() { return mpoles_; }
      // Returns bs1 of shell pair
      std::shared_ptr<BasisSet> bs1() { return bs1_; }
      // Returns bs2 of shell pair
      std::shared_ptr<BasisSet> bs2() { return bs2_; }
};

class PSI_API CFMMBox : public std::enable_shared_from_this<CFMMBox> {

    protected:
      // Parent of the CFMMBox
      std::weak_ptr<CFMMBox> parent_;
      // Children of the CFMMBox
      std::vector<std::shared_ptr<CFMMBox>> children_;

      // The primary shell pairs belonging to this box (empty if none)
      std::vector<std::shared_ptr<ShellPair>> primary_shell_pairs_;
      // The auxiliary shell pairs belonging to this box (empty if none)
      std::vector<std::shared_ptr<ShellPair>> auxiliary_shell_pairs_;

      // The box's origin (lower-left-front corner)
      Vector3 origin_;
      // Center of the box
      Vector3 center_;
      // Length of the box
      double length_;
      // Level the box is at (0 = root)
      int level_;
      // Maximum Multipole Angular Momentum
      int lmax_;
      // Well-separatedness criterion for this box
      int ws_;
      // Maximum well-separatedness for any given shell in the box 
      // (same as ws_ except for the most diffuse boxes in the level)
      int ws_max_;

      // Number of threads the calculation is running on
      int nthread_;

      // Multipoles of the box (Density-Matrix contracted), one for each density matrix (calculated for the ket basis)
      std::vector<std::shared_ptr<RealSolidHarmonics>> mpoles_;
      // Far field vector of the box, one for each density matrix (based on the multipoles of the ket basis)
      std::vector<std::shared_ptr<RealSolidHarmonics>> Vff_;

      // A list of all the near-field boxes to this box
      std::vector<std::shared_ptr<CFMMBox>> near_field_;
      // A list of all of the local-far-field boxes to this box
      std::vector<std::shared_ptr<CFMMBox>> local_far_field_;

      // Returns a shared pointer to the CFMMBox object
      std::shared_ptr<CFMMBox> get() { return shared_from_this(); }
      
    public:
      /// CFMMBox Constructor
      CFMMBox(std::shared_ptr<CFMMBox> parent, std::vector<std::shared_ptr<ShellPair>> primary_shell_pairs,
                std::vector<std::shared_ptr<ShellPair>> auxiliary_shell_pairs, Vector3 origin, double length,
                int level, int lmax, int ws);

      // Make children for this multipole box
      void make_children();
      // Sets the near field and local far field regions of the box
      void set_regions();

      // Compute multipoles for the box are contracted with the density matrix (depending on contraction type)
      void compute_multipoles(const std::vector<SharedMatrix>& D, ContractionType contraction_type);

      // Compute multipoles from children
      void compute_mpoles_from_children();
      // Computes the far field contribution from a far away sibling
      void compute_far_field_contribution(std::shared_ptr<CFMMBox> lff_box);
      // Compute the far field contibution from the parents
      void add_parent_far_field_contribution();

      // => USEFUL SETTER METHODS <= //

      // Set the maximum ws of the box
      void set_ws_max(int ws_max) { ws_max_ = ws_max; }

      // => USEFUL GETTER METHODS <= //

      // Get the multipole level the box is on
      int get_level() { return level_; }
      // Get the ws criterion of the box
      int get_ws() { return ws_; }
      // Get the value of a particular multipole (for the Nth density matrix)
      double get_mpole_val(int N, int l, int mu) { return mpoles_[N]->get_multipoles()[l][mu]; }
      // Get the far field value of a multipole (for the Nth density matrix)
      double get_Vff_val(int N, int l, int mu) { return Vff_[N]->get_multipoles()[l][mu]; }
      // Get the children of the box
      std::vector<std::shared_ptr<CFMMBox>>& get_children() { return children_; }
      // Get the bra shell pairs of the box
      std::vector<std::shared_ptr<ShellPair>>& get_primary_shell_pairs() { return primary_shell_pairs_; }
      // Get the ket shell pairs of the box
      std::vector<std::shared_ptr<ShellPair>>& get_auxiliary_shell_pairs() { return auxiliary_shell_pairs_; }
      // Gets the number of shell pairs in the box
      int nshell_pair() { return primary_shell_pairs_.size() + auxiliary_shell_pairs_.size(); }
      // Gets the number of bra shell pairs in the box
      int primary_nshell_pair() { return primary_shell_pairs_.size(); }
      // Gets the number of ket shell pairs in the box
      int auxiliary_nshell_pair() { return auxiliary_shell_pairs_.size(); }
      // Get the center of this box
      Vector3 center() { return center_; }
      // Gets the near_field_boxes of the box
      std::vector<std::shared_ptr<CFMMBox>>& near_field_boxes() { return near_field_; }
      // Gets the local far field boxes of the box
      std::vector<std::shared_ptr<CFMMBox>>& local_far_field_boxes() { return local_far_field_; }
      // Get the multipoles
      std::vector<std::shared_ptr<RealSolidHarmonics>>& multipoles() { return mpoles_; }
      // Gets the far field vector
      std::vector<std::shared_ptr<RealSolidHarmonics>>& far_field_vector() { return Vff_; }

}; // End class CFMMBox

class PSI_API CFMMTree {

    protected:
      // The molecule that this tree structure references
      std::shared_ptr<Molecule> molecule_;
      // The primary basis set
      std::shared_ptr<BasisSet> primary_;
      // The auxiliary basis set
      std::shared_ptr<BasisSet> auxiliary_;
      // List of all the significant primary shell-pairs in the molecule (U_SHELL, V_SHELL), U >= V
      std::vector<std::shared_ptr<ShellPair>> primary_shell_pairs_;
      // List of all the significant auxiliary shell-pairs in the molecule (SHELL, 0)
      std::vector<std::shared_ptr<ShellPair>> auxiliary_shell_pairs_;
      // What type of contraction is being performed? (Inferred by input parameters)
      ContractionType contraction_type_;

      // Number of Levels in the CFMM Tree
      int nlevels_;
      // Maximum Multipole Angular Momentum
      int lmax_;

      // The tree structure (implemented as list for random access)
      std::vector<std::shared_ptr<CFMMBox>> tree_;
      // List of all the leaf boxes (sorted by number of shell pairs for parallel efficiency)
      std::vector<std::shared_ptr<CFMMBox>> sorted_leaf_boxes_;
      // Harmonic Coefficients used to calculate multipoles
      std::shared_ptr<HarmonicCoefficients> mpole_coefs_;
      // Numerical cutoff for ERI screening
      double cutoff_;

      // Options object
      Options& options_;
      // Number of threads
      int nthread_;
      // Print flag, defaults to 1
      int print_;
      // Bench flag, defaults to 0
      int bench_;

      // List of all the primary shell-pairs to compute
      std::vector<std::pair<int, int>> primary_shellpair_tasks_;
      // Index from the primary shell-pair index to the bra shell pair
      std::vector<std::shared_ptr<ShellPair>> primary_shellpair_list_;
      // The box each primary shell-pair belongs to
      std::vector<std::shared_ptr<CFMMBox>> primary_shellpair_to_box_;
      // List of all the near field boxes that belong to a given primary shell-pair
      std::vector<std::vector<std::shared_ptr<CFMMBox>>> primary_shellpair_to_nf_boxes_;

      // List of all the auxiliary shell-pairs to compute
      std::vector<std::pair<int, int>> auxiliary_shellpair_tasks_;
      // Index from the auxiliary shell-pair index to the bra shell pair
      std::vector<std::shared_ptr<ShellPair>> auxiliary_shellpair_list_;
      // The box each auxiliary shell-pair belongs to
      std::vector<std::shared_ptr<CFMMBox>> auxiliary_shellpair_to_box_;
      // List of all the near field boxes that belong to a given auxiliary shell-pair
      std::vector<std::vector<std::shared_ptr<CFMMBox>>> auxiliary_shellpair_to_nf_boxes_;

      // local far-field box pairs at a given level of the tree
      std::vector<std::vector<std::pair<std::shared_ptr<CFMMBox>, std::shared_ptr<CFMMBox>>>> lff_task_pairs_per_level_;

      // Use density-based integral screening?
      bool density_screening_;
      // ERI Screening Tolerance
      double ints_tolerance_;

      // => Functions called ONLY once <= //

      // Make the root node of the CFMMTree
      void make_root_node();
      // Create children
      void make_children();
      // Sort the leaf nodes by number of shell-pairs
      void sort_leaf_boxes();
      // Set up near field and far field information for each box in the tree
      void setup_regions();
      // Setup shell-pair information and calculate multipoles for each shell-pair
      void setup_shellpair_info();
      // Set up information on local far field task pairs per level
      void setup_local_far_field_task_pairs();
      // Calculate the shell-pair multipoles at each leaf box (primary or auxiliary)
      void calculate_shellpair_multipoles(bool is_primary);

      // => Functions called ONCE per iteration <= //

      // Calculate multipoles
      void calculate_multipoles(const std::vector<SharedMatrix>& D);
      // Helper method to compute far field
      void compute_far_field();

      // Build near-field J (Gateway function, links to specific J builds based on contraction 
      void build_nf_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                      const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
		      const std::vector<double>& Jmet_max);
      // Build near-field J using Direct SCF algorithm (Jpq = (pq|rs)Drs)
      void build_nf_direct_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                      const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J);
      // Build gammaP's near field (gammaP = (P|uv)Duv)
      void build_nf_gamma_P(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                      const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
		      const std::vector<double>& Jmet_max);
      // Build density-fitted J's near field (Jpq = (pq|Q)*gammaQ)
      void build_nf_df_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                      const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
		      const std::vector<double>& Jmet_max);
      // Builds the near field interactions of the Coulomb metric with an auxiliary density
      void build_nf_metric(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                      const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J);

      // Build far-field J (long-range multipole interactions)
      void build_ff_J(std::vector<SharedMatrix>& J);

      // => ERI Screening <= //
      bool shell_significant(int P, int Q, int R, int S, std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                             const std::vector<SharedMatrix>& D);
    
    public:
      /// Constructor (automatically sets up the tree)
      /// Pass in null pointer if no primary or auxiliary
      CFMMTree(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);

      // Build the J matrix of CFMMTree
      void build_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, 
                    const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
		    const std::vector<double>& Jmet_max = {});
      // Returns the max tree depth
      int nlevels() { return nlevels_; }
      // Returns the max multipole AM
      int lmax() { return lmax_; }
      // Flip the contraction type (for DF integrals)
      void df_set_contraction(ContractionType contraction_type);
      // Print the CFMM Tree out
      void print_out();

}; // End class CFMMTree

} // namespace psi

#endif
