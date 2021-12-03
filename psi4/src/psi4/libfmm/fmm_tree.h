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

#define ERFCI10 (4.572824967389485)

namespace psi {

class Options;

class PSI_API ShellPair {
    protected:
      // The basisset associated with the shell-pair
      std::shared_ptr<BasisSet> basisset_;
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
      ShellPair(std::shared_ptr<BasisSet>& basisset, std::pair<int, int> pair_index, std::shared_ptr<HarmonicCoefficients>& mpole_coefs);

      // Calculate the multipole moments of the Shell-Pair about a center
      void calculate_mpoles(Vector3 box_center, std::shared_ptr<OneBodyAOInt> s_ints,
                            std::shared_ptr<OneBodyAOInt> mpole_ints, int lmax);

      // Returns the shell pair index
      std::pair<int, int> get_shell_pair_index() { return pair_index_; }
      // Returns the center of the shell pair
      Vector3 get_center() { return center_; }
      // Returns the radial extent of the shell pair
      double get_extent() { return extent_; }
      // Returns the multipole moments of the shell pairs about a center
      std::vector<std::shared_ptr<RealSolidHarmonics>>& get_mpoles() { return mpoles_; }
};

class PSI_API CFMMBox : public std::enable_shared_from_this<CFMMBox> {

    protected:
      // Parent of the CFMMBox
      std::weak_ptr<CFMMBox> parent_;
      // Children of the CFMMBox
      std::vector<std::shared_ptr<CFMMBox>> children_;

      // The shell pairs belonging to this box
      std::vector<std::shared_ptr<ShellPair>> shell_pairs_;

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

      // Number of threads the calculation is running on
      int nthread_;

      // Multipoles of the box (Density-Matrix contracted)
      std::shared_ptr<RealSolidHarmonics> mpoles_;
      // Far field vector of the box
      std::shared_ptr<RealSolidHarmonics> Vff_;

      // A list of all the near-field boxes to this box
      std::vector<std::shared_ptr<CFMMBox>> near_field_;
      // A list of all of the local-far-field boxes to this box
      std::vector<std::shared_ptr<CFMMBox>> local_far_field_;

      // Returns a shared pointer to the CFMMBox object
      std::shared_ptr<CFMMBox> get() { return shared_from_this(); }
      
    public:
      // Generic Constructor
      CFMMBox(std::shared_ptr<CFMMBox> parent, std::vector<std::shared_ptr<ShellPair>> shell_pairs, 
              Vector3 origin, double length, int level, int lmax, int ws);

      // Make children for this multipole box
      void make_children();
      // Compute multipoles directly
      void compute_mpoles(std::shared_ptr<BasisSet>& basisset, std::vector<SharedMatrix>& D);
      // Compute multipoles from children
      void compute_mpoles_from_children();
      // Sets the near field and local far field and calculates far field vector from local and parent far fields
      void compute_far_field();

      // => USEFUL GETTER METHODS <= //
      // Get the multipole level the box is on
      int get_level() { return level_; }
      // Get the ws criterion of the box
      int get_ws() { return ws_; }
      // Get the value of a particular multipole
      double get_mpole_val(int l, int mu) { return mpoles_->get_multipoles()[l][mu]; }
      // Get the far field value of a multipole
      double get_Vff_val(int l, int mu) { return Vff_->get_multipoles()[l][mu]; }
      // Get the children of the box
      std::vector<std::shared_ptr<CFMMBox>>& get_children() { return children_; }
      // Get the shell pairs of the box
      std::vector<std::shared_ptr<ShellPair>>& get_shell_pairs() { return shell_pairs_; }
      // Gets the number of shell pairs in the box
      int get_nsp() { return shell_pairs_.size(); }
      // Gets the near_field_boxes of the box
      std::vector<std::shared_ptr<CFMMBox>>& near_field_boxes() { return near_field_; }
      // Gets the far field vector
      std::shared_ptr<RealSolidHarmonics>& far_field_vector() { return Vff_; }

}; // End class CFMMBox

class PSI_API CFMMTree {

    protected:
      // The molecule that this tree structure references
      std::shared_ptr<Molecule> molecule_;
      // The basis set that the molecule uses
      std::shared_ptr<BasisSet> basisset_;
      // Density Matrix of Molecule
      std::vector<SharedMatrix> D_;
      // Coulomb Matrix of Molecule
      std::vector<SharedMatrix> J_;
      // List of all the significant shell-pairs in the molecule
      std::vector<std::shared_ptr<ShellPair>> shell_pairs_;
      // Number of Levels in the CFMM Tree
      int nlevels_;
      // Maximum Multipole Angular Momentum
      int lmax_;
      // The tree structure (implemented as list for random access)
      std::vector<std::shared_ptr<CFMMBox>> tree_;
      // Harmonic Coefficients used to calculate multipoles
      std::shared_ptr<HarmonicCoefficients> mpole_coefs_;
      // A list of significant bra shell pairs
      std::vector<std::pair<int, int>> nf_bra_shell_pairs_;

      // Options object
      Options& options_;
      // Number of threads
      int nthread_;

      // The integral objects used to compute the integrals
      std::vector<std::shared_ptr<TwoBodyAOInt>> ints_;

      // Sort the shell-pairs (radix sort)
      void sort_shell_pairs();
      // Make the root node of the CFMMTree
      void make_root_node();
      // Create children
      void make_children();
      // Calculate multipoles
      void calculate_multipoles();
      // Helper method to compute far field
      void compute_far_field();
      // Build near-field J (like Direct SCF)
      void build_nf_J();
      // Build far-field J (long-range multipole interactions)
      void build_ff_J();
    
    public:
      // Constructor
      CFMMTree(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, std::vector<SharedMatrix>& D, 
                std::vector<SharedMatrix>& J, Options& options);

      // Build the J matrix of CFMMTree
      void build_J();
      // Print the CFMM Tree out
      void print_out();

}; // End class CFMMTree

} // namespace psi

#endif
