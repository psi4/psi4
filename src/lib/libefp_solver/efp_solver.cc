/**
 * EFP solver
 */ 

#include <boost/regex.hpp>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility> 
             
#include <libmints/mints.h>
#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
//#include <liboptions/python.h>
#include <psifiles.h>
#include <libmints/multipolesymmetry.h>
#include <psi4-dec.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "efp_solver.h"
#include "../../libefp/libefp/src/efp.h"

using namespace boost;
using namespace std;
using namespace psi;

regex efpAtomSymbol("A\\d*([A-Z]{1,2})\\d*", regbase::normal | regbase::icase);
smatch reMatches;


// TODO: change allocated memory to shared pointers and ditch the deletes

namespace psi { namespace efp {
EFP::EFP(Options& options): options_(options) 
{
fprintf(outfile,"efp constructor calling common_init\n");
	common_init();
}

EFP::~EFP(){
	efp_shutdown(efp_);
}

static bool is_lib(const char *name)
{
	size_t len = strlen(name);

	return name[len - 2] == '_' && (name[len - 1] == 'l' || name[len - 1] == 'L');
}

static size_t name_len(const char *name)
{
	return is_lib(name) ? strlen(name) - 2 : strlen(name);
}

/**
 * Basic creation of EFP object and reading of options
 */
void EFP::common_init() {
    enum efp_result res;
    int print = options_.get_int("PRINT");

//fprintf(outfile,"efp::common_init reads efp options from local options_\n");
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    //elst_enabled_ = options_.get_bool("EFP_ELST");
    //pol_enabled_  = options_.get_bool("EFP_POL");
    //disp_enabled_ = options_.get_bool("EFP_DISP");
    //exch_enabled_ = options_.get_bool("EFP_EXCH");

    //std::string dertype = options_.get_str("DERTYPE");

//fprintf(outfile, "\nDERTYPE: %s\n", dertype.c_str());
    //do_grad_ = false;
    //if (dertype == "FIRST")
    //    do_grad_ = true;

    //if (elst_enabled_)
    //    opts.terms |= EFP_TERM_ELEC;
    //if (pol_enabled_)
    //    opts.terms |= EFP_TERM_POL;
    //if (disp_enabled_)
    //    opts.terms |= EFP_TERM_DISP;
    //if (exch_enabled_)
    //    opts.terms |= EFP_TERM_XR;

    //std::string elst_damping = options_.get_str("EFP_ELST_DAMPING");
    //std::string disp_damping = options_.get_str("EFP_DISP_DAMPING");

    //if (elst_damping == "SCREEN")
    //    opts.elec_damp = EFP_ELEC_DAMP_SCREEN;
    //else if (elst_damping == "OVERLAP")
    //    opts.elec_damp = EFP_ELEC_DAMP_OVERLAP;
    //else if (elst_damping == "OFF")
    //    opts.elec_damp = EFP_ELEC_DAMP_OFF;

    //if (disp_damping == "TT")
    //    opts.disp_damp = EFP_DISP_DAMP_TT;
    //else if (disp_damping == "OVERLAP")
    //    opts.disp_damp = EFP_DISP_DAMP_OVERLAP;
    //else if (disp_damping == "OFF")
    //    opts.disp_damp = EFP_DISP_DAMP_OFF;

    //molecule_ = Process::environment.molecule();
    //nfrag_ = options_["FRAGS"].size();

//fprintf(outfile,"efp::common_init calls efp_create\n");
    efp_ = efp_create();

    if (!efp_)
        throw PsiException("EFP::common_init():", __FILE__, __LINE__);

//fprintf(outfile,"efp::common_init calls efp_set_opts\n");
    if (res = efp_set_opts(efp_, &opts))
        throw PsiException("EFP::common_init(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    //fprintf(outfile, "\n\n");
    //fprintf(outfile, "%s", efp_banner());
    //fprintf(outfile, "\n\n");
    //fprintf(outfile, "  ==> Calculation Information <==\n\n");
    //fprintf(outfile, "  Electrostatics damping: %12s\n", elst_damping.c_str());

    //if (disp_enabled_)
    //    fprintf(outfile, "  Dispersion damping:     %12s\n", disp_damping.c_str());

    //fprintf(outfile, "  Electrostatics enabled:  %12d\n", elst_enabled_);
    //fprintf(outfile, "  Polarization enabled:    %12d\n", pol_enabled_);
    //fprintf(outfile, "  Dispersion enabled:      %12d\n", disp_enabled_);
    //fprintf(outfile, "  Exchanged enabled:       %12d\n", exch_enabled_);
    //fprintf(outfile, "  Gradient enabled:        %12d\n", do_grad_);
    //fprintf(outfile, "\n");

    //fprintf(outfile, "\n");
}

/**
 * Add potential files and names for all fragments
 */
void EFP::add_fragments(std::vector<std::string> fnames)
{
    enum efp_result res;
	int n_uniq;
	char path[256];
	std::vector<std::string> uniq;

    std::string psi_data_dir = Process::environment("PSIDATADIR");
    std::string fraglib_path = psi_data_dir + "/fraglib";
    std::string userlib_path = ".";

    nfrag_ = fnames.size();

    for (int i = 0; i < fnames.size(); i++) {
        std::string name = fnames[i];
		std::transform(name.begin(), name.end(), name.begin(), ::tolower);
		uniq.push_back(name);
	}

	std::sort(uniq.begin(), uniq.end());
	n_uniq = 1;

    for (int i = 1; i < fnames.size(); i++)
		if (uniq[i - 1] != uniq[i])
			uniq[n_uniq++] = uniq[i];

	for (int i = 0; i < n_uniq; i++) {
		std::string name = uniq[i];
		const char *prefix = is_lib(name.c_str()) ? fraglib_path.c_str() : userlib_path.c_str();

		strcat(strncat(strcat(strcpy(path, prefix), "/"), name.c_str(),
			name_len(name.c_str())), ".efp");
		if (res = efp_add_potential(efp_, path))
            throw PsiException("EFP::add_fragments(): " + std::string (efp_result_to_string(res)) + name,__FILE__,__LINE__);
	}

    for (int i = 0; i < fnames.size(); i++)
        if (res = efp_add_fragment(efp_, fnames[i].c_str()))
            throw PsiException("EFP::add_fragments(): " + std::string (efp_result_to_string(res)) + " " + fnames[i],__FILE__,__LINE__);
}

/**
 * Get gradient of the interaction energy of the EFP electrostatics with the QM nuclei (point charges)
 */
boost::shared_ptr<Vector> EFP::get_electrostatic_gradient() {

    int natom = molecule_->natom();
    boost::shared_ptr<Vector> grad ( new Vector( 3*natom ) );
    double * grad_p = grad->pointer();

    // verify natom matches the number of point charges in efp
    enum efp_result res;
    int n_ptc;
    if ( res = efp_get_point_charge_count(efp_,&n_ptc) ) {
        throw PsiException("EFP::get_electrostatic_gradient(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    }
    if ( n_ptc != natom ) {
        throw PsiException("EFP::get_electrostatic_gradient(): natom does not match number of point charges in efp",__FILE__,__LINE__);
    }
    
    if ( res = efp_get_point_charge_gradient(efp_,grad_p) ) {
        throw PsiException("EFP::get_electrostatic_gradient(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    }

    return grad;
}

/**
 * Provide list of coordinates of quantum mechanical atoms to efp_set_point_charges
 */
void EFP::set_qm_atoms(){

    int natom = molecule_->natom();

    boost::shared_ptr<Vector>   q (new Vector(natom));
    boost::shared_ptr<Vector> xyz (new Vector(3*natom));

    double * q_p   = q->pointer();
    double * xyz_p = xyz->pointer();
    for (int A = 0; A < natom; A++) {
        if ( molecule_->Z(A) == 0.0 ) continue;
        q_p[A]       = molecule_->Z(A);
        xyz_p[3*A]   = molecule_->x(A);
        xyz_p[3*A+1] = molecule_->y(A);
        xyz_p[3*A+2] = molecule_->z(A);
    }

    enum efp_result res;
    if ( res = efp_set_point_charges(efp_,natom,q_p,xyz_p) ) {
        throw PsiException("EFP::SetQMAtoms(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    }
}

/**
 * Add fragment by name
 */
void EFP::add_fragment(std::string fname) {
    enum efp_result res;

    if (res = efp_add_fragment(efp_, fname.c_str()))
        throw PsiException("EFP::add_fragment(): " + std::string (efp_result_to_string(res)) + " " + fname,__FILE__,__LINE__);
}

/**
 * Set points or xyzabc coordinates for all fragments simultaneously
 */
void EFP::set_coordinates(int type, double * coords) {
    enum efp_result res;
    enum efp_coord_type ctype;

    if(type == 0)
        ctype = EFP_COORD_TYPE_XYZABC;
    else if(type == 1)
        ctype = EFP_COORD_TYPE_POINTS;
    else if(type == 2)
        ctype = EFP_COORD_TYPE_ROTMAT;

    if ((res = efp_set_coordinates(efp_, ctype, coords)))
        throw PsiException("EFP::set_coordinates(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

/*
 * Set points or xyzabc coordinates for all atoms in a fragment
 */
void EFP::set_frag_coordinates(int frag_idx, int type, double * coords) {
    enum efp_result res;
    enum efp_coord_type ctype;

    if(type == 0)
        ctype = EFP_COORD_TYPE_XYZABC;
    else if(type == 1)
        ctype = EFP_COORD_TYPE_POINTS;
    else if(type == 2)
        ctype = EFP_COORD_TYPE_ROTMAT;

    int local;
    efp_get_frag_count(efp_, &local);

    if ((res = efp_set_frag_coordinates(efp_, frag_idx, ctype, coords)))
        throw PsiException("EFP::set_frag_coordinates(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

// this function returns a shared matrix containing the efp contribution to the potential
// felt by qm atoms in an scf procedure.
boost::shared_ptr<Matrix> EFP::modify_Fock() {

    // get number of multipoles (charges, dipoles, quadrupoles, octupoles)
    int n_multipole = 0;
    if ( efp_get_multipole_count(efp_,&n_multipole) != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return number of multipoles",__FILE__,__LINE__);
    }

    // multipole coordinates are stored array xyz.
    boost::shared_ptr<Vector> xyz  (new Vector(3*n_multipole));
    boost::shared_ptr<Vector> mult (new Vector((1+3+6+10)*n_multipole));

    // get multipoles from libefp
    // dipoles stored as     x,y,z
    // quadrupoles stored as xx,yy,zz,xy,xz,yz
    // octupoles stored as   xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
    if ( efp_get_multipole_coordinates(efp_,xyz->pointer()) != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return multipole coordinates",__FILE__,__LINE__);
    }
    if ( efp_get_multipole_values(efp_,mult->pointer()) != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return multipole values",__FILE__,__LINE__);
    }

    // induced dipoles
    int n_id = 0;
    if ( efp_get_induced_dipole_count(efp_,&n_id)  != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return number of induced dipoles",__FILE__,__LINE__);
    }
    boost::shared_ptr<Vector> xyz_id (new Vector(3*n_id));
    if ( efp_get_induced_dipole_coordinates(efp_,xyz_id->pointer()) != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return induced dipole coordinates",__FILE__,__LINE__);
    }
    boost::shared_ptr<Vector> id (new Vector(3*n_id));
    if ( efp_get_induced_dipole_values(efp_,id->pointer()) != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return induced dipole values",__FILE__,__LINE__);
    }
    boost::shared_ptr<Vector> idt (new Vector(3*n_id));
    if ( efp_get_induced_dipole_conj_values(efp_,idt->pointer())  != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return induced dipole conjugate values",__FILE__,__LINE__);
    }
    // take average of induced dipole and conjugate
    id->add(idt);
    id->scale(0.5);

//    // normal multipole integrals are ordered as follows 
//    // x, y, z, 
//    // xx, xy, xz, yy, yz, zz, 
//    // xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
//    // presumably the new integrals will be ordered similarly
//
//    // arrays to map our multipole ordering to Ilya's
//    int mapq[6]  = { 0, 3, 4, 1, 5, 2};
//    int mapo[10] = { 0, 3, 4, 5, 9, 7, 1, 6, 8, 2};

    // contract/dot/something multipoles with multipole integrals.  the result goes into V
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<OneBodyAOInt> efp_ints(wfn->integral()->ao_efp_multipole_potential());

    int nbf = wfn->basisset()->nbf();

                               // 0    X    Y    Z      XX       YY       ZZ       XY       XZ       YZ
    const double prefacs[20] = { 1.0, 1.0, 1.0, 1.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0,
    //   XXX       YYY       ZZZ       XXY       XXZ       XYY       YYZ       XZZ       YZZ       XYZ  
      1.0/15.0, 1.0/15.0, 1.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 6.0/15.0};

    std::vector<SharedMatrix> mats;
    for(int i=0; i < 20; ++i) {
        mats.push_back(SharedMatrix(new Matrix(nbf, nbf)));
    }

    SharedMatrix V(new Matrix("EFP contribution to the Fock Matrix", nbf, nbf));

    // e-Nu contributions from fragments ( this is now done below with the other multipoles ).
    //V->add( EFP_nuclear_potential() );

    // multipole contributions to Fock matrix
    double * xyz_p  = xyz->pointer();
    double * mult_p = mult->pointer();
    for (int n = 0; n < n_multipole; n++) {
        for(int i=0; i < 20; ++i){
           mats[i]->zero();
        }
        Vector3 coords(xyz_p[n*3],xyz_p[n*3+1],xyz_p[n*3+2]);
        efp_ints->set_origin(coords);
        efp_ints->compute(mats);

        // add point charges from atoms to multipoles at atom center
        for (int frag = 0; frag < nfrag_; frag++) {
            int natom = 0;
            if ( efp_get_frag_atom_count(efp_,frag,&natom) != EFP_RESULT_SUCCESS ) {
                throw PsiException("libefp failed to return the number of atoms",__FILE__,__LINE__);
            }
            efp_atom * atoms = (efp_atom*)malloc(natom*sizeof(efp_atom));
            if ( efp_get_frag_atoms(efp_, frag, natom, atoms) != EFP_RESULT_SUCCESS ) {
                throw PsiException("libefp failed to return atom charges",__FILE__,__LINE__);
            }
            for (int i = 0; i < natom; i++) {
                double dx = atoms[i].x - xyz_p[n*3];
                if ( fabs(dx) > 1e-10 ) continue;
                double dy = atoms[i].y - xyz_p[n*3+1];
                if ( fabs(dy) > 1e-10 ) continue;
                double dz = atoms[i].z - xyz_p[n*3+2];
                if ( fabs(dz) > 1e-10 ) continue;

                mult_p[20*n] += atoms[i].znuc;
                break;
            }
        }

        for(int i=0; i < 20; ++i){
            mats[i]->scale( -prefacs[i] * mult_p[20*n+i] );
            V->add(mats[i]);
        }
    }

    // induced dipole contributions to Fock matrix
    xyz_p  = xyz_id->pointer();
    mult_p = id->pointer();
    for (int n = 0; n < n_id; n++) {
        for(int i=0; i < 20; ++i){
           mats[i]->zero();
        }
        Vector3 coords(xyz_p[n*3],xyz_p[n*3+1],xyz_p[n*3+2]);
        efp_ints->set_origin(coords);
        efp_ints->compute(mats);
        // only dealing with dipoles here:
        for(int i=0; i < 3; ++i){
            mats[i+1]->scale( -prefacs[i+1] * mult_p[3*n+i] );
            V->add(mats[i+1]);
        }
    }

    V->print();
    return V;
}

double EFP::EFP_QM_nuclear_repulsion_energy() {
    double nu = 0.0;
    boost::shared_ptr<Molecule> mol = Process::environment.molecule();
    for (int frag = 0; frag < nfrag_; frag++) {
        int natom = 0;
        if ( efp_get_frag_atom_count(efp_,frag,&natom) != EFP_RESULT_SUCCESS ) {
            throw PsiException("libefp failed to return the number of atoms",__FILE__,__LINE__);
        }
        efp_atom * atoms = (efp_atom*)malloc(natom*sizeof(efp_atom));
        if ( efp_get_frag_atoms(efp_, frag, natom, atoms) != EFP_RESULT_SUCCESS ) {
            throw PsiException("libefp failed to return atom charges",__FILE__,__LINE__);
        }

        for (int i = 0; i < natom; i++) {
            double znuc = atoms[i].znuc;
            double x    = atoms[i].x;
            double y    = atoms[i].y;
            double z    = atoms[i].z;

            for (int j = 0; j < mol->natom(); j++) {
                double dx = x - mol->x(j);
                double dy = y - mol->y(j);
                double dz = z - mol->z(j);
                double r  = sqrt(dx*dx+dy*dy+dz*dz);
                nu += znuc * mol->Z(j) / r;
            }


        }
    }
    return nu;
}

/**
 * compute efp contribution to scf energy
 */
double EFP::scf_energy_update() {
    enum efp_result res;
    double efp_energy;

    if (res = efp_get_wavefunction_dependent_energy(efp_, &efp_energy))
        throw PsiException("EFP::scf_energy_update(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return efp_energy;
}

/**
 * Compute and return efp total energy only
 */
double EFP::ComputeEnergy() {
    enum efp_result res;

    // Main EFP computation routine 
    if (res = efp_compute(efp_, do_grad_ ? 1 : 0))
        throw PsiException("EFP::Compute():efp_compute() " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    
    struct efp_energy energy;
    
    if (res = efp_get_energy(efp_, &energy))
        throw PsiException("EFP::Compute():efp_get_energy(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return energy.total;
}
/**
 * Compute efp energy components and/or gradient
 */
void EFP::Compute() {
    enum efp_result res;

    // Main EFP computation routine 
    if (res = efp_compute(efp_, do_grad_ ? 1 : 0))
        throw PsiException("EFP::Compute():efp_compute() " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    
    struct efp_energy energy;
    
    if (res = efp_get_energy(efp_, &energy))
        throw PsiException("EFP::Compute():efp_get_energy(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    
    if (do_grad_) {
        SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag_, 6));
        double ** psmgrad = smgrad->pointer();
        if (res = efp_get_gradient(efp_, nfrag_, psmgrad[0]))
            throw PsiException("EFP::Compute():efp_get_gradient(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

        fprintf(outfile, "  ==> EFP Gradient <==\n\n");

        for (int i = 0; i < nfrag_; i++) {
            for (int j = 0; j < 6; j++) {
                fprintf(outfile, "%14.6lf", psmgrad[i][j]);
            }
            fprintf(outfile, "\n");
        }
        fprintf(outfile, "\n");

//        torque_ = smgrad;
//        psi::Process::environment.set_efp_torque(smgrad);
        psi::Process::environment.set_gradient(smgrad);
        //smgrad->print();
    }

    fprintf(outfile, "  ==> Energetics <==\n\n");

    fprintf(outfile, "  EFP/EFP Electrostatics Energy = %20.12f [H] %s\n", 
        energy.electrostatic + energy.charge_penetration, elst_enabled_ ? "*" : "");
    if (do_qm_) {
        fprintf(outfile, "  QM-Nuc/EFP Electrostatics Energy = %20.12f [H] %s\n", 
            energy.electrostatic_point_charges, qm_elst_enabled_ ? "*" : "");
        fprintf(outfile, "  Total Electrostatics Energy = %20.12f [H] %s\n", 
            energy.electrostatic + energy.charge_penetration + energy.electrostatic_point_charges, 
            (elst_enabled_ || qm_elst_enabled_) ? "*" : "");
    }
    fprintf(outfile, "  %7s Polarization Energy =   %20.12f [H] %s\n", do_qm_ ? "" : "EFP/EFP", 
        energy.polarization, pol_enabled_ ? "*" : "");
    fprintf(outfile, "  EFP/EFP Dispersion Energy =     %20.12f [H] %s\n", 
        energy.dispersion, disp_enabled_ ? "*" : "");
    fprintf(outfile, "  EFP/EFP Exchange Energy =       %20.12f [H] %s\n", 
        energy.exchange_repulsion, exch_enabled_ ? "*" : "");
    fprintf(outfile, "  Total Energy =                  %20.12f [H] %s\n", 
        energy.total, "*");

    Process::environment.globals["EFP ELST ENERGY"] = energy.electrostatic + energy.charge_penetration + energy.electrostatic_point_charges;
    Process::environment.globals["EFP POL ENERGY"] = energy.polarization;
    Process::environment.globals["EFP DISP ENERGY"] = energy.dispersion;
    Process::environment.globals["EFP EXCH ENERGY"] = energy.exchange_repulsion;
    Process::environment.globals["CURRENT ENERGY"] = energy.total;
}

/**
 * Get number of fragments
 */
int EFP::get_frag_count(void) {
    enum efp_result res;
    int n=0;

    if (res = efp_get_frag_count(efp_, &n))
        throw PsiException("EFP::get_frag_count(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    nfrag_ = n;
    return n;
}

/**
 * Get number of atoms in a fragment
 */
int EFP::get_frag_atom_count(int frag_idx) {
    enum efp_result res;
    int n=0;

    if (res = efp_get_frag_atom_count(efp_, frag_idx, &n))
        throw PsiException("EFP::get_frag_atom_count(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return n;
}

/**
 * Get charge on a fragment
 */
double EFP::get_frag_charge(int frag_idx) {
    enum efp_result res;
    double charge=0.0;

    if (res = efp_get_frag_charge(efp_, frag_idx, &charge))
        throw PsiException("EFP::get_frag_charge(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return charge;
}

/**
 * Get multiplicity for a fragment
 */
int EFP::get_frag_multiplicity(int frag_idx) {
    enum efp_result res;
    int multiplicity=0;

    if (res = efp_get_frag_multiplicity(efp_, frag_idx, &multiplicity))
        throw PsiException("EFP::get_frag_multiplicity(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return multiplicity;
}

/**
 * Get atomic numbers of atoms in a fragment
 * NOTE: not sure if we need all these input checks here
 *   as long as don't plan on making invalid calls,
 *   can trap libefp errors instead and throw error
 *   rather than returning NULL. Same for other fns below.
 */
double *EFP::get_frag_atom_Z(int frag_idx) {

    if (frag_idx >= nfrag_) return NULL;

    int frag_natom = get_frag_atom_count(frag_idx);
    if (frag_natom == 0) return NULL;

    struct efp_atom atoms[frag_natom];
    enum efp_result res;
    res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms);

    if (res != EFP_RESULT_SUCCESS) return NULL;

    double *frag_atom_Z = new double[frag_natom];
    for (int i=0; i<frag_natom; ++i)
        frag_atom_Z[i] = atoms[i].znuc;

    return frag_atom_Z;
}

/**
 * Get the mass of all atoms in a fragment
 */
double *EFP::get_frag_atom_mass(int frag_idx) {
  
    if (frag_idx >= nfrag_) return NULL;
  
    int frag_natom = get_frag_atom_count(frag_idx);
  
    struct efp_atom atoms[frag_natom];
    enum efp_result res;
    res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms);
    if (res != EFP_RESULT_SUCCESS) return NULL;
  
    double *frag_atom_mass = new double[frag_natom];
    for (int i=0; i<frag_natom; ++i)
        frag_atom_mass[i] = atoms[i].mass;
  
    return frag_atom_mass;
}

/**
 * Get COM for a fragment from libefp coordinates
 */
double *EFP::get_com(int frag_idx) {
    if (frag_idx >= nfrag_) return NULL;
  
    double *xyzabc = new double [6*nfrag_];
    efp_get_coordinates(efp_, nfrag_, xyzabc);

    double *com = new double[3];
    com[0] = xyzabc[6*frag_idx+0];
    com[1] = xyzabc[6*frag_idx+1];
    com[2] = xyzabc[6*frag_idx+2];
  
    return com;
}

/**
 * Get the xyz coordinates of all atoms in a fragment
 */
double *EFP::get_frag_atom_coord(int frag_idx) {

    if (frag_idx >= nfrag_) return NULL;
    
    int frag_natom = get_frag_atom_count(frag_idx);
    if (frag_natom == 0) return NULL;
    
    struct efp_atom atoms[frag_natom];
    enum efp_result res;
    res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms);
    
    if (res != EFP_RESULT_SUCCESS) return NULL;
    
    double *frag_atom_coord = new double[3*frag_natom];
    for (int i=0; i<frag_natom; ++i) {
        frag_atom_coord[3*i]   = atoms[i].x;
        frag_atom_coord[3*i+1] = atoms[i].y;
        frag_atom_coord[3*i+2] = atoms[i].z;
    }

    return frag_atom_coord;
}

/**
 * Get the atom label of all atoms in a fragment
 */
std::vector<std::string> EFP::get_frag_atom_label(int frag_idx) {

    if (frag_idx >= nfrag_)
        throw PsiException("EFP::get_frag_atom_label():",__FILE__,__LINE__);

    int frag_natom = get_frag_atom_count(frag_idx);
    if (frag_natom == 0)
        throw PsiException("EFP::get_frag_atom_label():",__FILE__,__LINE__);

    struct efp_atom atoms[frag_natom];
    enum efp_result res;
    res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms);

    if (res != EFP_RESULT_SUCCESS)
        throw PsiException("EFP::get_frag_atom_label():",__FILE__,__LINE__);

    std::vector<std::string> frag_atom_label;
    for (int i=0; i<frag_natom; ++i) {
        std::string label = atoms[i].label;
        frag_atom_label.push_back(label);
    }

    return frag_atom_label;
}

/**
 * Print simple private members of class
 */
void EFP::print_out() {

    fprintf(outfile, "  ==> EFP/EFP Setup <==\n\n");
    fprintf(outfile, "  Number of EFP fragments: %12d\n", nfrag_);
    fprintf(outfile, "  Electrostatics enabled?: %12s\n", elst_enabled_ ? "true" : "false");
    fprintf(outfile, "  Polarization enabled?:   %12s\n", pol_enabled_ ? "true" : "false");
    fprintf(outfile, "  Dispersion enabled?:     %12s\n", disp_enabled_ ? "true" : "false");
    fprintf(outfile, "  Exchange enabled?:       %12s\n", exch_enabled_ ? "true" : "false");
    if (elst_enabled_)
        fprintf(outfile, "  Electrostatics damping:  %12s\n", elst_damping_.c_str());
    if (pol_enabled_)
        fprintf(outfile, "  Polarization damping:    %12s\n", pol_damping_.c_str());
    if (disp_enabled_)
        fprintf(outfile, "  Dispersion damping:      %12s\n", disp_damping_.c_str());
    fprintf(outfile, "  Gradient enabled?:       %12s\n", do_grad_ ? "true" : "false");

    if (do_qm_) {
        fprintf(outfile, "\n  ==> QM/EFP Setup <==\n\n");
        fprintf(outfile, "  Number of QM fragments:  %12d\n", molecule_->nfragments());
        fprintf(outfile, "  Electrostatics enabled?: %12s\n", qm_elst_enabled_ ? "true" : "false");
        fprintf(outfile, "  Polarization enabled?:   %12s\n", qm_pol_enabled_ ? "true" : "false");
        fprintf(outfile, "  Dispersion enabled?:     %12s\n", "undefined");
        fprintf(outfile, "  Exchange enabled?:       %12s\n", "undefined");
    }

    print_efp_geometry();

    if (do_qm_) {
        fprintf(outfile,"  ==> QM Geometry <==\n\n");
        molecule_->print();
    }
}


efp_result electron_density_field_fn(int n_pt, const double *xyz, double *field, void *user_data) {
    // These should all be members of the SCF class in the final implementation.
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Molecule> mol = wfn->molecule();
    boost::shared_ptr<BasisSet> basis = wfn->basisset();
    boost::shared_ptr<OneBodyAOInt> field_ints(wfn->integral()->electric_field());

    int nbf = basis->nbf();
    std::vector<SharedMatrix> intmats;
    intmats.push_back(SharedMatrix(new Matrix("Ex integrals", nbf, nbf)));
    intmats.push_back(SharedMatrix(new Matrix("Ey integrals", nbf, nbf)));
    intmats.push_back(SharedMatrix(new Matrix("Ez integrals", nbf, nbf)));

    SharedMatrix Da = wfn->Da_subset("AO");
    SharedMatrix Db;
    if(!wfn->same_a_b_orbs())
        Db = wfn->Db_subset("AO");

    for(int n = 0; n < n_pt; ++n){
        field_ints->set_origin(Vector3(xyz[3*n], xyz[3*n+1],xyz[3*n+2]));
        for(int m = 0; m < 3; ++m)
            intmats[m]->zero();
        field_ints->compute(intmats);
        double Ex = Da->vector_dot(intmats[0]);
        double Ey = Da->vector_dot(intmats[1]);
        double Ez = Da->vector_dot(intmats[2]);
        if(wfn->same_a_b_dens()){
            Ex *= 2.0;
            Ey *= 2.0;
            Ez *= 2.0;
        }else{
            Ex += Db->vector_dot(intmats[0]);
            Ey += Db->vector_dot(intmats[1]);
            Ez += Db->vector_dot(intmats[2]);
        }
        // Vector3 nucterms = ElectricFieldInt::nuclear_contribution(field_ints->origin(), mol);
        // Ex += nucterms[0];
        // Ey += nucterms[1];
        // Ez += nucterms[2];
        field[3*n]   = Ex;
        field[3*n+1] = Ey;
        field[3*n+2] = Ez;
    }
    return EFP_RESULT_SUCCESS;
}

/**
 * Resetting of EFP options
 */
void EFP::set_options() {

    molecule_ = Process::environment.molecule();

    //fprintf(outfile,"efp::set_options() passing options to libefp\n");

    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    elst_enabled_    = options_.get_bool("EFP_ELST");
    pol_enabled_     = options_.get_bool("EFP_POL");
    disp_enabled_    = options_.get_bool("EFP_DISP");
    exch_enabled_    = options_.get_bool("EFP_EXCH");
    do_qm_           = options_.get_bool("QMEFP");
    qm_elst_enabled_ = do_qm_ && options_.get_bool("QMEFP_ELST");
    qm_pol_enabled_  = do_qm_ && options_.get_bool("QMEFP_POL");

    std::string dertype = options_.get_str("DERTYPE");
    do_grad_ = false;
    if (dertype == "FIRST")
        do_grad_ = true;

    // AI_DISP, AI_XR, AI_CHTR may be enabled in a future libefp release
    if (qm_elst_enabled_)
        opts.terms |= EFP_TERM_AI_ELEC;
    if (qm_pol_enabled_)
        opts.terms |= EFP_TERM_AI_POL;

    // CHTR may be enabled in a future libefp release
    if (elst_enabled_)
        opts.terms |= EFP_TERM_ELEC;
    if (pol_enabled_)
        opts.terms |= EFP_TERM_POL;
    if (disp_enabled_)
        opts.terms |= EFP_TERM_DISP;
    if (exch_enabled_)
        opts.terms |= EFP_TERM_XR;

    elst_damping_ = options_.get_str("EFP_ELST_DAMPING");
    disp_damping_ = options_.get_str("EFP_DISP_DAMPING");
    pol_damping_ = options_.get_str("EFP_POL_DAMPING");

    if (elst_damping_ == "SCREEN")
        opts.elec_damp = EFP_ELEC_DAMP_SCREEN;
    else if (elst_damping_ == "OVERLAP")
        opts.elec_damp = EFP_ELEC_DAMP_OVERLAP;
    else if (elst_damping_ == "OFF")
        opts.elec_damp = EFP_ELEC_DAMP_OFF;

    if (disp_damping_ == "TT")
        opts.disp_damp = EFP_DISP_DAMP_TT;
    else if (disp_damping_ == "OVERLAP")
        opts.disp_damp = EFP_DISP_DAMP_OVERLAP;
    else if (disp_damping_ == "OFF")
        opts.disp_damp = EFP_DISP_DAMP_OFF;

    if (pol_damping_ == "TT")
        opts.pol_damp = EFP_POL_DAMP_TT;
    else if (pol_damping_ == "OFF")
        opts.pol_damp = EFP_POL_DAMP_OFF;

    enum efp_result res;
    if (res = efp_set_opts(efp_, &opts))
        throw PsiException("EFP::set_options(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    fprintf(outfile, "\n\n");
    fprintf(outfile, "%s", efp_banner());
    fprintf(outfile, "\n\n");
    //fprintf(outfile, "  ==> Calculation Information (This section going away) <==\n\n");
    //fprintf(outfile, "  Electrostatics damping: %12s\n", elst_damping_.c_str());

    //if (disp_enabled_)
    //    fprintf(outfile, "  Dispersion damping:     %12s\n", disp_damping_.c_str());

    //fprintf(outfile, "  Electrostatics enabled:  %12d\n", elst_enabled_);
    //fprintf(outfile, "  Polarization enabled:    %12d\n", pol_enabled_);
    //fprintf(outfile, "  Dispersion enabled:      %12d\n", disp_enabled_);
    //fprintf(outfile, "  Exchange enabled:        %12d\n", exch_enabled_);
    //fprintf(outfile, "  Gradient enabled:        %12d\n", do_grad_);
    //fprintf(outfile, "\n");

    //fprintf(outfile, "\n");

    // sets call-back function to provide electric field from electrons
    efp_set_electron_density_field_fn( efp_, electron_density_field_fn );

}

int EFP::efp_natom() {
    if ( nfrag_ == 0 ) return 0;
    int natom = 0;
    for (int frag = 0; frag < nfrag_; frag++) {
        std::vector<std::string> symbol = get_frag_atom_label(frag);
        natom += symbol.size();
    }
    return natom;
}

void EFP::print_efp_geometry() {
    if ( nfrag_ == 0 ) return;

    Molecule::GeometryUnits units = Process::environment.molecule()->units();
    bool is_angstrom = ( units == Molecule::Angstrom ) ;

    fprintf(outfile,"\n");
    fprintf(outfile,"  ==> EFP Geometry <==\n\n");
    fprintf(outfile,"    Geometry (in %s):\n\n", //, charge = %d, multiplicity = %d:\n\n",
            is_angstrom ? "Angstrom" : "Bohr"); //, get_frag_multiplicity, multiplicity_);
    fprintf(outfile,"     Center            X                  Y                  Z        \n");
    fprintf(outfile,"    --------   -----------------  -----------------  -----------------\n");

    for (int frag = 0; frag < nfrag_; frag++) {
        double * xyz = get_frag_atom_coord(frag);
        std::vector<std::string> symbol = get_frag_atom_label(frag);
        int natom = symbol.size();
        for (int i = 0; i < natom; i++) {
            regex_match(symbol[i], reMatches, efpAtomSymbol);
            fprintf(outfile,"   %5s     ",reMatches[1].str().c_str(),frag);
            fflush(outfile);
            for (int j = 0; j < 3; j++) {
                fprintf(outfile,"  %17.12lf", xyz[i*3+j] * ( is_angstrom ? pc_bohr2angstroms : 1.0 ) );
            }
            fprintf(outfile," (EFP %3i)\n",frag+1);
            fflush(outfile);
        }

    }
    fprintf(outfile,"\n");
}

boost::shared_ptr<Matrix> EFP::EFP_nuclear_potential() {

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

    int nbf = wfn->basisset()->nbf();

    boost::shared_ptr<Matrix> V (new Matrix(nbf,nbf));

    for (int frag = 0; frag < nfrag_; frag++) {
        int natom = 0;
        if ( efp_get_frag_atom_count(efp_,frag,&natom) != EFP_RESULT_SUCCESS ) {
            throw PsiException("libefp failed to return the number of atoms",__FILE__,__LINE__);
        }
        efp_atom * atoms = (efp_atom*)malloc(natom*sizeof(efp_atom));
        if ( efp_get_frag_atoms(efp_, frag, natom, atoms) != EFP_RESULT_SUCCESS ) {
            throw PsiException("libefp failed to return atom charges",__FILE__,__LINE__);
        }

        SharedMatrix V_charge(new Matrix("External Potential (Charges)", nbf, nbf));

        SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", natom, 4));
        double** Zxyzp = Zxyz->pointer();
        for (int i = 0; i < natom; i++) {
            Zxyzp[i][0] = atoms[i].znuc;
            Zxyzp[i][1] = atoms[i].x;
            Zxyzp[i][2] = atoms[i].y;
            Zxyzp[i][3] = atoms[i].z;
        }

        boost::shared_ptr<PotentialInt> pot(static_cast<PotentialInt*>(wfn->integral()->ao_potential()));
        pot->set_charge_field(Zxyz);
        pot->compute(V_charge);

        // add to efp fock contribution
        V->add(V_charge);

        V_charge.reset();
        pot.reset();
    }
    return V;
}


}

} // End namespaces

