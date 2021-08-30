#ifndef libfmm_fmm_tree_H
#define libfmm_fmm_tree_H

#include "psi4/pragma.h"

#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libfmm/multipoles_helper.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <unordered_map>

#define ERFCI10 (4.572824967389485)

namespace psi {

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
      std::shared_ptr<CFMMBox> get();

      // Compute the near field and far field J matrix contributions
      void compute_nf_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J);
      void compute_ff_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J);
      
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
      // Sets the near field and local far field vectors
      void set_nf_lff();
      // Calculate far field vector from local and parent far fields
      void compute_far_field_vector();
      // Compute the box's contribution to the J matrix
      void compute_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J);

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
      std::vector<std::shared_ptr<ShellPair>>& get_shell_pairs() { return shell_pairs_; };

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
      // The tree structure (implemented as list for simplification)
      std::vector<std::shared_ptr<CFMMBox>> tree_;
      // Harmonic Coefficients used to calculate multipoles
      std::shared_ptr<HarmonicCoefficients> mpole_coefs_;

      // Sort the shell-pairs (radix sort)
      void sort_shell_pairs();
      // Make the root node of the CFMMTree
      void make_root_node();
      // Create children
      void make_children();
      // Calculate multipoles
      void calculate_multipoles();
      // Helper method to set the near field and lff vectors
      void set_nf_lff();
      // Helper method to compute far field
      void compute_far_field();
      // Helper method to build the J Matrix recursively
      void calculate_J(CFMMBox* box);
    
    public:
      // Constructor
      CFMMTree(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, 
                std::vector<SharedMatrix>& J, const std::vector<std::pair<int, int>>& shell_pairs, int nlevels, int lmax);

      // Build the J matrix of CFMMTree
      void build_J();
      // Print the CFMM Tree out
      void print_out();

}; // End class CFMMTree

} // namespace psi

#endif