#ifndef EFP_SOLVER_H
#define EFP_SOLVER_H
/*
 * EFP header
 */

#include<libmints/molecule.h>
#include "../../interfaces/libefp/libefp/src/efp.h"

struct efp;

namespace psi {
    class Options;
}

namespace boost {
    template<class T> class shared_ptr;
}

namespace psi { 

class Vector;

namespace efp {


class EFP {
    // warning: options_ is pointer to current options object and may not reflect
    //    proper efp options outside of common_init()
    Options & options_;
    protected:

        /// Number of fragments
        int nfrag_;

        /// EFP struct for libefp
        struct efp * efp_;

        /// active QM molecule
        boost::shared_ptr<Molecule> molecule_;

        /// Flags for EFP/EFP options
        bool elst_enabled_, pol_enabled_, disp_enabled_, exch_enabled_;
        std::string elst_damping_, pol_damping_, disp_damping_;

        /// Flags for QM/EFP options
        bool qm_elst_enabled_, qm_pol_enabled_;

        /// Flags for flow control options
        bool do_grad_, do_qm_;

        /// Initialize options
        void common_init();

        /// If a gradient is available it will be here:
        SharedMatrix torque_;

    public:
        /// Constructor
        EFP(Options& options);

        /// Destructor
        ~EFP();

        /// Add potential files and names for all fragments
        void add_fragments(std::vector<std::string> fnames);

        /// Sets the geometry hints for a given fragment; wrapper to efp_set_frag_coordinates
        void set_frag_coordinates(int frag_idx, int type, double * coords);

        /// Finalize fragment composition of efp system; wrapper to efp_prepare
        void finalize_fragments();

        /// Returns multiplicity for a given fragment; wrapper to efp_get_frag_multiplicity
        int get_frag_multiplicity(int frag_idx);

        /// Returns charge for a given fragment; wrapper to efp_get_frag_charge
        double get_frag_charge(int frag_idx);

        /// Returns the number of atoms in specified fragment; wrapper to efp_get_frag_atom_count
        int get_frag_atom_count(int frag_idx);

        /// Returns atomic numbers of all atoms in a given fragment
        double * get_frag_atom_Z(int frag_idx);

        /// Returns masses of all atoms in a given fragment
        double * get_frag_atom_mass(int frag_idx);

        /// Returns atom label of all atoms in a given fragment
        std::vector<std::string> get_frag_atom_label(int frag_idx);

        /// Returns xyz coordinates of all atoms in a given fragment
        double * get_frag_atom_coord(int frag_idx);

        /// Sets member data and libefp options to Options values
        ///    Usually called through py_psi_efp_set_options so correct module Options prepared
        void set_options();

        /// Designate which atoms are QM; wrapper to efp_set_point_charges
        void set_qm_atoms();

        /// Returns EFP permanent moment contribution to V
        boost::shared_ptr<Matrix> modify_Fock_permanent();

        /// Returns EFP induced dipole contribution to V
        boost::shared_ptr<Matrix> modify_Fock_induced();

        /// Returns EFP contribution to SCF energy; wrapper to efp_get_wavefunction_dependent_energy
        double scf_energy_update();

        /// Compute energy and/or gradient
        void compute();

        /// Returns the number of EFP fragments; wrapper to efp_get_frag_count
        int get_frag_count(void);

        /// Print all of the EFP atoms
        void print_efp_geometry();

        /// Prints private members of efp object
        void print_out(void);

        /// Sets the EFP gradient
        void set_torque(const SharedMatrix& torq) { torque_ = torq; }

        /// Returns the gradient
        SharedMatrix torque() const { return torque_; }
};

}
}

//        /// Returns the center of mass of a given fragment
//        double * get_com(int frag_idx);
//        /// Number of EFP atoms
//        int efp_natom();
//        /// Computes the nuclear potential integrals from the EFP fragments
//        boost::shared_ptr<Matrix> EFP_nuclear_potential();
//        /// Wrapper to efp_get_point_charge_gradient
//        boost::shared_ptr<Vector> get_electrostatic_gradient();
//        /// Add EFP fragment
//        void add_fragment(std::string fname);
//        /// Sets the geometry hints for all fragments at once
//        void set_coordinates(int type, double * coords);
//        /// Computes the nuclear repulsion between the QM and EFP regions
//        double EFP_QM_nuclear_repulsion_energy();

#endif
