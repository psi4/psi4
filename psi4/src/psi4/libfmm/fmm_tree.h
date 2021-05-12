#ifndef libfmm_fmm_tree_H
#define libfmm_fmm_tree_H

#include "psi4/pragma.h"

#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libfmm/multipoles_helper.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <unordered_map>

#define ERFCI10 (4.572824967389485)

namespace psi {

class CFMMBox {

    protected:
      // Parent of the CFMMBox
      CFMMBox* parent_;
      // Children of the CFMMBox
      std::vector<CFMMBox *> children_;
      // Level the box is at (0 = root)
      int level_;
      // Maximum Multipole Angular Momentum
      int lmax_;
      // Well-separatedness criteria for the particular Box
      int ws_;
      // The molecule that is referenced by this box
      std::shared_ptr<Molecule> molecule_;
      // The basis set that the molecule uses
      std::shared_ptr<BasisSet> basisset_;
      // Density Matrix of Molecule
      std::vector<SharedMatrix> D_;
      // A reference to the Coulomb Matrix of the molecule (every box can modify it)
      std::vector<SharedMatrix> J_;
      // The atoms in the molecule which are centered within the bounds of the box
      std::vector<int> atoms_;
      // Length of the box
      double length_;
      // The box's origin (lower-left-front corner)
      Vector3 origin_;
      // Center of the box
      Vector3 center_;

      // Number of threads the calculation is running on
      int nthread_;

      // Solid Harmonics Coefficients of Box
      std::shared_ptr<HarmonicCoefficients> mpole_coefs_;

      // Multipoles of the box, Density Contracted
      std::shared_ptr<RealSolidHarmonics> mpoles_;
      // Far field vector of the box, Density Contracted
      std::shared_ptr<RealSolidHarmonics> Vff_;

      // A list of all the near-field boxes to this box
      std::vector<CFMMBox *> near_field_;
      // A list of all of the local-far-field boxes to this box
      std::vector<CFMMBox *> local_far_field_;
      // Far field energy
      double ff_energy_;

      // Common function used by constructor
      void common_init(CFMMBox* parent, std::shared_ptr<Molecule> molecule, 
                        std::shared_ptr<BasisSet> basisset, Vector3 origin, double length, int level, int lmax);

      // Compute the J matrix contributions at each level
      void compute_self_J();
      void compute_nf_J();
      void compute_ff_J();
      
    public:
      // Constructor for a root box
      CFMMBox(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, 
                std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J, int lmax);
      // Constructor for child boxes
      CFMMBox(CFMMBox* parent, std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, 
                std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J, Vector3 origin, double length, int level, int lmax);
      // Make children for this multipole box
      void make_children();
      // Compute multipoles directly
      void compute_mpoles();
      // Compute multipoles from children
      void compute_mpoles_from_children();
      // Sets the near field and local far field vectors
      void set_nf_lff();
      // Calculate far field vector from local and parent far fields
      void compute_far_field_vector();
      // Compute the box's contribution to the J matrix
      void compute_J();

      // => USEFUL GETTER METHODS <= //

      // Get the multipole level the box is on
      int get_level() { return level_; }
      // Get the number of atoms in the box
      int natom() { return atoms_.size(); }
      // Get the value of a particular multipole
      double get_mpole_val(int l, int mu) { return mpoles_->get_multipoles()[l][mu]; }
      // Get the far field value of a multipole
      double get_Vff_val(int l, int mu) { return Vff_->get_multipoles()[l][mu]; }
      // Get the children of the box
      std::vector<CFMMBox*>& get_children() { return children_; }
      // Return the far field energy
      double ff_energy() { return ff_energy_; }

      // Destructor
      virtual ~CFMMBox();
      

}; // End class CFMMBox

class CFMMTree {

    protected:
      // Number of Levels in the CFMM Tree
      int nlevels_;
      // The molecule that this tree structure references
      std::shared_ptr<Molecule> molecule_;
      // The basis set that the molecule uses
      std::shared_ptr<BasisSet> basisset_;
      // Maximum Multipole Angular Momentum
      int lmax_;
      // Root of this tree structure
      CFMMBox* root_;
      // Density Matrix of Molecule
      std::vector<SharedMatrix> D_;
      // Coulomb Matrix of Molecule
      std::vector<SharedMatrix> J_;
      // Far field energy
      double ff_energy_;

      // Create children
      void make_children(CFMMBox* box);
      // Calculate multipoles
      void calculate_multipoles(CFMMBox* box);
      // Helper method to set the near field and lff vectors
      void set_nf_lff(CFMMBox* box);
      // Helper method to compute far field
      void compute_far_field(CFMMBox* box);
      // Helper method to build the J Matrix recursively
      void calculate_J(CFMMBox* box);
    
    public:
      // Constructor
      CFMMTree(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, 
                std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J, int nlevels, int lmax);

      // Build the J matrix of CFMMTree
      void build_J();

      // Destructor
      virtual ~CFMMTree();

}; // End class CFMMTree

} // namespace psi

#endif