#ifndef EFP_SOLVER_H
#define EFP_SOLVER_H
/*
 * EFP header
 */

#include<libmints/molecule.h>

struct efp;

namespace psi{
  class Options;
}

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace efp{


class EFP {
    // warning: options_ is pointer to current options object, and may not reflect
    // proper efp options outside of common_init()
    Options & options_;
    protected:
        int nfrag_;
        struct efp * efp_;
        boost::shared_ptr<Molecule>molecule_;
        bool elst_enabled_, pol_enabled_, disp_enabled_, exch_enabled_, do_grad_;
        /// Initialize options
        void common_init();
    public:
        EFP(Options& options);
        ~EFP();
  
        /// Set geometry
        void SetGeometry();

        /// Compute energy and/or gradietn
        void Compute();

        /// Designate which atoms are qm atoms
        void SetQMAtoms();

        /// Returns EFP contribution to SCF energy
        double scf_energy_update();

        /// Returns EFP contribution to V
        boost::shared_ptr<Matrix> modify_Fock();

	    // Make list of frag names
	    char *make_name_list();

	    // Make list of .efp file names
	    char *make_potential_file_list(const char *, const char *, const char *);

        /// Returns the number of EFP fragments
        int get_frag_count(void);
        
        /// Returns the number of atoms in specified fragment
        int get_frag_atom_count(int frag_idx);
        
        /// Returns atomic numbers of all atoms in a given fragment
        double *get_frag_atom_Z(int frag_idx);
        
        /// Returns masses of all atoms in a given fragment
        double *get_frag_atom_mass(int frag_idx);

        /// Returns the center of mass of a given fragment
        double *get_com(int frag_idx);

        ///
        void set_nfragments(int nfrag) { nfrag_ = nfrag; }

        ///
        int get_nfragments(void) { return nfrag_; }

        /// Sets the geometry hints for all fragments at once
        void set_coordinates(int type, double * coords);

        /// Sets the geometry hints for a given fragment
        void set_frag_coordinates(int frag_idx, int type, double * coords);

        /// Returns xyz coordinates of all atoms in a given fragment
        double *get_frag_atom_coord(int frag_idx);

        /// Returns atom label of all atoms in a given fragment
        std::vector<std::string> get_frag_atom_label(int frag_idx);
};

}
}

#endif
