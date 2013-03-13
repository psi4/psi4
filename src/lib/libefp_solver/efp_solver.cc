/**
 * EFP solver
 */ 

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
#include <liboptions/python.h>
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

// TODO: change allocated memory to shared pointers and ditch the deletes

namespace psi { namespace efp {
EFP::EFP(Options& options): options_(options) 
{
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

    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    elst_enabled_ = options_.get_bool("EFP_ELST");
    pol_enabled_  = options_.get_bool("EFP_POL");
    disp_enabled_ = options_.get_bool("EFP_DISP");
    exch_enabled_ = options_.get_bool("EFP_EXCH");

    std::string dertype = options_.get_str("DERTYPE");
    do_grad_ = false;
    if (dertype == "FIRST")
        do_grad_ = true;

    if (elst_enabled_)
        opts.terms |= EFP_TERM_ELEC;
    if (pol_enabled_)
        opts.terms |= EFP_TERM_POL;
    if (disp_enabled_)
        opts.terms |= EFP_TERM_DISP;
    if (exch_enabled_)
        opts.terms |= EFP_TERM_XR;

    std::string elst_damping = options_.get_str("EFP_ELST_DAMPING");
    std::string disp_damping = options_.get_str("EFP_DISP_DAMPING");

    if (elst_damping == "SCREEN")
        opts.elec_damp = EFP_ELEC_DAMP_SCREEN;
    else if (elst_damping == "OVERLAP")
        opts.elec_damp = EFP_ELEC_DAMP_OVERLAP;
    else if (elst_damping == "OFF")
        opts.elec_damp = EFP_ELEC_DAMP_OFF;

    if (disp_damping == "TT")
        opts.disp_damp = EFP_DISP_DAMP_TT;
    else if (disp_damping == "OVERLAP")
        opts.disp_damp = EFP_DISP_DAMP_OVERLAP;
    else if (disp_damping == "OFF")
        opts.disp_damp = EFP_DISP_DAMP_OFF;

    molecule_ = Process::environment.molecule();

    efp_ = efp_create();

    if (!efp_)
        throw PsiException("EFP::common_init():", __FILE__, __LINE__);

    if (res = efp_set_opts(efp_, &opts))
        throw PsiException("EFP::common_init(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    fprintf(outfile, "\n\n");
    fprintf(outfile, "%s", efp_banner());
    fprintf(outfile, "\n\n");
    fprintf(outfile, "  ==> Calculation Information <==\n\n");
    fprintf(outfile, "  Electrostatics damping: %12s\n", elst_damping.c_str());

    if (disp_enabled_)
        fprintf(outfile, "  Dispersion damping:     %12s\n", disp_damping.c_str());

    fprintf(outfile, "  Electrostatics enabled:  %12d\n", elst_enabled_);
    fprintf(outfile, "  Polarization enabled:    %12d\n", pol_enabled_);
    fprintf(outfile, "  Dispersion enabled:      %12d\n", disp_enabled_);
    fprintf(outfile, "  Exchanged enabled:       %12d\n", exch_enabled_);
    fprintf(outfile, "  Gradient enabled:        %12d\n", do_grad_);
    fprintf(outfile, "\n");

    fprintf(outfile, "\n");
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
 * Provide list of coordinates of quantum mechanical atoms
 */
void EFP::SetQMAtoms(){
    //TODO: extend molecule class and coordentry class to separate qm and efp atoms
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

    // get number of multipoles
    int * n_multipole = (int*)malloc(4*sizeof(int));
    enum efp_result err = efp_get_multipole_count(efp_,n_multipole);
    if ( err != EFP_RESULT_SUCCESS ) {
        throw PsiException("libefp failed to return number of multipoles",__FILE__,__LINE__);
    }

    // workspace for efp_get_multipoles.
    //double ** xyz = (double**)malloc(4*sizeof(double*));
    //double **   z = (double**)malloc(4*sizeof(double*));
    //for (int i = 0; i < 4; i++) {
    //    xyz[i] = (double*)malloc(3*n_multipole[i]*sizeof(double));
    //    for (int j = 0; j < 3*n_multipole[i]; j++) xyz[i][j] = 0.0;
    //}
    //z[0] = (double*)malloc(n_multipole[0]*sizeof(double));
    //z[1] = (double*)malloc(3*n_multipole[1]*sizeof(double));
    //z[2] = (double*)malloc(6*n_multipole[2]*sizeof(double));
    //z[3] = (double*)malloc(10*n_multipole[3]*sizeof(double));

    //for (int j = 0; j < n_multipole[0]; j++)    z[0][j] = 0.0;
    //for (int j = 0; j < 3*n_multipole[1]; j++)  z[1][j] = 0.0;
    //for (int j = 0; j < 6*n_multipole[2]; j++)  z[2][j] = 0.0;
    //for (int j = 0; j < 10*n_multipole[3]; j++) z[3][j] = 0.0;

    // get multipoles from libefp
    // dipoles stored as     x,y,z
    // quadrupoles stored as xx,yy,zz,xy,xz,yz
    // octupoles stored as   xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
    //
    // err = efp_get_grag_atoms - for atom charges
    // err = efp_get_multipole_values - for electrostatics multipones
    // err = efp_get_induced_dipole_values - for polarization induced dipoles
    // err = efp_get_induced_dipole_conj_values - for polarization induced dipoles
    //
    // induced dipoles
    //
    // int n_id;
    // check_fail(efp_get_induced_dipole_count(impl_->efp, &n_id));
    // double *xyz_id = new double[n_id * 3];
    // check_fail(efp_get_induced_dipole_coordinates(impl_->efp, xyz_id));
    // double *id = new double[n_id * 3];
    // check_fail(efp_get_induced_dipole_values(impl_->efp, id));
    // double *idt = new double[n_id * 3];
    // check_fail(efp_get_induced_dipole_conj_values(impl_->efp, idt));
    //
    // take avarage of id and idt, 0.5*(id+idt)

    // get electrostatic potential at each point returned in the xyz array
    // TODO: need this function

    // grab matrix factory from wavefunction
    boost::shared_ptr<Wavefunction> wfn         = Process::environment.wavefunction();
    boost::shared_ptr<MatrixFactory> matrix     = wfn->matrix_factory();
    boost::shared_ptr<IntegralFactory> integral = wfn->integral();

    // generate multipole integrals:
    // 
    // they will be ordered as follows in the vector, multipoles
    // x, y, z, 
    // xx, xy, xz, yy, yz, zz, 
    // xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
    // 
    boost::shared_ptr<OneBodySOInt> mult3 ( integral->so_multipoles(3) );
    boost::shared_ptr<MultipoleSymmetry> multsym ( new MultipoleSymmetry(3,molecule_,integral,matrix) );
    std::vector<boost::shared_ptr<Matrix> > multipoles = multsym->create_matrices("Multipole: ");
    mult3->compute(multipoles);

    // arrays to map our multipole ordering to Ilya's
    int mapq[6]  = { 0, 3, 4, 1, 5, 2};
    int mapo[10] = { 0, 3, 4, 5, 9, 7, 1, 6, 8, 2};

    // dot multipoles with multipole integrals.  the result goes into V
    // TODO: need this function
    boost::shared_ptr<Matrix> V = matrix->create_shared_matrix("EFP V contribution");

    // free workspace memory needed by libefp
    //free(n_multipole);
    //for (int i = 0; i < 4; i++) {
    //    free(xyz[i]);
    //    free(z[i]);
    //}
    //free(xyz);
    //free(z);
   
    return V;
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
 * Compute efp energy components and/or gradient
 */
void EFP::Compute() {
    enum efp_result res;
    double *grad = NULL;

    // Main EFP computation routine 
    if (res = efp_compute(efp_, do_grad_ ? 1 : 0))
        throw PsiException("EFP::Compute():efp_compute() " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    
    struct efp_energy energy;
    
    if (res = efp_get_energy(efp_, &energy))
        throw PsiException("EFP::Compute():efp_get_energy(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    
    if (do_grad_) {
        grad = new double[6 * nfrag_];
        if (res = efp_get_gradient(efp_, nfrag_, grad))
            throw PsiException("EFP::Compute():efp_get_gradient(): " + std::string (efp_result_to_string(res)),__FILE__,__LINE__);

        fprintf(outfile, "  ==> EFP Gradient <==\n\n");

        double *pgrad = grad;
        for (int i = 0; i < nfrag_; i++) {
            for (int j = 0; j < 6; j++) {
                fprintf(outfile, "%14.6lf", *pgrad++);
            }
            fprintf(outfile, "\n");
        }
        fprintf(outfile, "\n");

        SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag_, 6));
        double ** psmgrad = smgrad->pointer();
        pgrad = grad;
        for (int i = 0; i < nfrag_; i++) {
            for (int jj = 0; jj < 6; jj++) {
                psmgrad[i][jj] = *pgrad++;
            }
        }

        psi::Process::environment.set_gradient(smgrad);
        smgrad->print();
    }

    fprintf(outfile, "  ==> Energetics <==\n\n");

    fprintf(outfile, "  Electrostatics Energy = %24.16f [H] %s\n", energy.electrostatic, elst_enabled_ ? "*" : "");
    fprintf(outfile, "  Polarization Energy =   %24.16f [H] %s\n", energy.polarization, pol_enabled_ ? "*" : "");
    fprintf(outfile, "  Dispersion Energy =     %24.16f [H] %s\n", energy.dispersion, disp_enabled_ ? "*" : "");
    fprintf(outfile, "  Exchange Energy =       %24.16f [H] %s\n", energy.exchange_repulsion, exch_enabled_ ? "*" : "");
    fprintf(outfile, "  Total Energy =          %24.16f [H] %s\n", energy.total, "*");

    Process::environment.globals["EFP ELST ENERGY"] = energy.electrostatic;
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

    fprintf(outfile, "  Number of efp fragments: %12d\n", nfrag_);
    fprintf(outfile, "  Electrostatics enabled:  %12d\n", elst_enabled_);
    fprintf(outfile, "  Polarization enabled:    %12d\n", pol_enabled_);
    fprintf(outfile, "  Dispersion enabled:      %12d\n", disp_enabled_);
    fprintf(outfile, "  Exchanged enabled:       %12d\n", exch_enabled_);
    fprintf(outfile, "  Gradient enabled:        %12d\n", do_grad_);
    fprintf(outfile, "\n");

    //molecule_->print();  // TODO: used to work, broken now?
}

}

} // End namespaces

