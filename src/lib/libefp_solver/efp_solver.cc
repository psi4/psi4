/*
 * EFP solver
 */

#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

#include <libmints/mints.h>
#include <physconst.h>
#include <psi4-dec.h>

#include "efp_solver.h"

boost::regex efpAtomSymbol("A\\d*([A-Z]{1,2})\\d*", boost::regbase::normal | boost::regbase::icase);
boost::smatch reMatches;

// TODO: change allocated memory to shared pointers and ditch the deletes


namespace psi { namespace efp {

EFP::EFP(Options& options): options_(options)
{
    common_init();
}

EFP::~EFP()
{
    efp_shutdown(efp_);
}

/*
 * Basic creation of EFP object and options structure
 */
void EFP::common_init()
{
    enum efp_result res;
    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    if (!(efp_ = efp_create()))
        throw PsiException("EFP::common_init():efp_create()", __FILE__, __LINE__);

    if ((res = efp_set_opts(efp_, &opts)))
        throw PsiException("EFP::common_init():efp_set_opts() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

/*
 * Add potential files and names for all fragments
 */
void EFP::add_fragments(std::vector<std::string> fnames)
{
    enum efp_result res;
    bool not_found = true;
    std::vector<std::string> uniq;

    // Paths to search for efp files: here + PSIPATH + library
    std::string libraryPath = Process::environment("PSIDATADIR") + "/fraglib";
    std::string efpPath = ".:" + Process::environment("PSIPATH") + ":" + libraryPath;
    boost::char_separator<char> sep(":");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(efpPath, sep);

    nfrag_ = fnames.size();

    for (size_t i=0; i<fnames.size(); i++) {
        std::string name = fnames[i];
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        uniq.push_back(name);
    }

    std::sort(uniq.begin(), uniq.end());
    int n_uniq = 1;

    for (size_t i=1; i<fnames.size(); i++)
        if (uniq[i - 1] != uniq[i])
            uniq[n_uniq++] = uniq[i];

    // Loop over unique fragments
    for (int i=0; i<n_uniq; i++) {
        std::string name = uniq[i];
        not_found = true;

        // Loop over possible locations
        for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            std::string psiPathWithFragment = *tok_iter + "/" + name + ".efp";
            boost::filesystem::path bf_path = boost::filesystem::system_complete(psiPathWithFragment);

            if (!(res = efp_add_potential(efp_, bf_path.string().c_str()))) {
                outfile->Printf("  EFP fragment %s read from %s\n",
                    name.c_str(), bf_path.string().c_str());
                not_found = false;
                break;
            }
        }
        if (not_found)
            throw PsiException("EFP::add_fragments(): Fragment " + name +
                " not located in " + efpPath +
                "\nLast error is " + std::string (efp_result_to_string(res)),
                __FILE__,__LINE__);
    }

    // Initialize each fragment (not just unique)
    for (size_t i=0; i<fnames.size(); i++)
        if ((res = efp_add_fragment(efp_, fnames[i].c_str())))
            throw PsiException("EFP::add_fragments(): " +
                std::string (efp_result_to_string(res)) + " " + fnames[i],__FILE__,__LINE__);
}

/*
 * Set points or xyzabc coordinates for all atoms in a fragment
 */
void EFP::set_frag_coordinates(int frag_idx, int type, double * coords)
{
    enum efp_result res;
    enum efp_coord_type ctype;

    if(type == 0)
        ctype = EFP_COORD_TYPE_XYZABC;
    else if(type == 1)
        ctype = EFP_COORD_TYPE_POINTS;
    else if(type == 2)
        ctype = EFP_COORD_TYPE_ROTMAT;

    size_t local;
    efp_get_frag_count(efp_, &local);

    if ((res = efp_set_frag_coordinates(efp_, frag_idx, ctype, coords)))
        throw PsiException("EFP::set_frag_coordinates() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

/*
 * Prepare libefp for calculation after all fragments added
 */
void EFP::finalize_fragments()
{
    enum efp_result res;

    if ((res = efp_prepare(efp_)))
        throw PsiException("EFP::finalize_fragments() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

/*
 * Get multiplicity for a fragment
 */
int EFP::get_frag_multiplicity(int frag_idx)
{
    enum efp_result res;
    int multiplicity=0;

    if ((res = efp_get_frag_multiplicity(efp_, frag_idx, &multiplicity)))
        throw PsiException("EFP::get_frag_multiplicity(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return multiplicity;
}

/*
 * Get charge on a fragment
 */
double EFP::get_frag_charge(int frag_idx)
{
    enum efp_result res;
    double charge=0.0;

    if ((res = efp_get_frag_charge(efp_, frag_idx, &charge)))
        throw PsiException("EFP::get_frag_charge(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return charge;
}

/*
 * Get number of atoms in a fragment
 */
int EFP::get_frag_atom_count(int frag_idx)
{
    enum efp_result res;
    size_t n=0;

    if ((res = efp_get_frag_atom_count(efp_, frag_idx, &n)))
        throw PsiException("EFP::get_frag_atom_count(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return n;
}

/*
 * Get atomic numbers of atoms in a fragment
 */
double * EFP::get_frag_atom_Z(int frag_idx)
{
    enum efp_result res;

    size_t frag_natom = get_frag_atom_count(frag_idx);
    struct efp_atom atoms[frag_natom];
    if ((res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms)))
        throw PsiException("EFP::get_frag_atom_Z() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    double *frag_atom_Z = new double[frag_natom];
    for (size_t i=0; i<frag_natom; ++i)
        frag_atom_Z[i] = atoms[i].znuc;

    return frag_atom_Z;
}

/*
 * Get the mass of all atoms in a fragment
 */
double * EFP::get_frag_atom_mass(int frag_idx)
{
    enum efp_result res;

    size_t frag_natom = get_frag_atom_count(frag_idx);
    struct efp_atom atoms[frag_natom];
    if ((res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms)))
        throw PsiException("EFP::get_frag_atom_mass() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    double *frag_atom_mass = new double[frag_natom];
    for (size_t i=0; i<frag_natom; ++i)
        frag_atom_mass[i] = atoms[i].mass;

    return frag_atom_mass;
}

/*
 * Get the atom label of all atoms in a fragment
 */
std::vector<std::string> EFP::get_frag_atom_label(int frag_idx)
{
    enum efp_result res;

    size_t frag_natom = get_frag_atom_count(frag_idx);
    struct efp_atom atoms[frag_natom];
    if ((res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms)))
        throw PsiException("EFP::get_frag_atom_label() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    std::vector<std::string> frag_atom_label;
    for (size_t i=0; i<frag_natom; ++i) {
        std::string label = atoms[i].label;
        frag_atom_label.push_back(label);
    }

    return frag_atom_label;
}

/*
 * Get the xyz coordinates of all atoms in a fragment
 */
double * EFP::get_frag_atom_coord(int frag_idx)
{
    enum efp_result res;

    size_t frag_natom = get_frag_atom_count(frag_idx);
    struct efp_atom atoms[frag_natom];
    if ((res = efp_get_frag_atoms(efp_, frag_idx, frag_natom, atoms)))
        throw PsiException("EFP::get_frag_atom_coord() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    double *frag_atom_coord = new double[3*frag_natom];
    for (size_t i=0; i<frag_natom; ++i) {
        frag_atom_coord[3*i]   = atoms[i].x;
        frag_atom_coord[3*i+1] = atoms[i].y;
        frag_atom_coord[3*i+2] = atoms[i].z;
    }

    return frag_atom_coord;
}

efp_result electron_density_field_fn(size_t n_pt, const double *xyz, double *field, void *user_data)
{
    // TODO These should all be members of the SCF class in the final implementation.
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Molecule> mol = wfn->molecule();
    boost::shared_ptr<BasisSet> basis = wfn->basisset();
    boost::shared_ptr<OneBodyAOInt> field_ints(wfn->integral()->electric_field());

    int nao = basis->nao();
    std::vector<SharedMatrix> intmats;
    intmats.push_back(SharedMatrix(new Matrix("Ex integrals", nao, nao)));
    intmats.push_back(SharedMatrix(new Matrix("Ey integrals", nao, nao)));
    intmats.push_back(SharedMatrix(new Matrix("Ez integrals", nao, nao)));

    SharedMatrix Da = wfn->Da_subset("AO");
    SharedMatrix Db;
    if (!wfn->same_a_b_orbs())
        Db = wfn->Db_subset("AO");

    for (size_t n=0; n<n_pt; ++n) {
        field_ints->set_origin(Vector3(xyz[3*n], xyz[3*n+1], xyz[3*n+2]));
        for (int m=0; m<3; ++m)
            intmats[m]->zero();
        field_ints->compute(intmats);
        double Ex = Da->vector_dot(intmats[0]);
        double Ey = Da->vector_dot(intmats[1]);
        double Ez = Da->vector_dot(intmats[2]);
        if (wfn->same_a_b_dens()) {
            Ex *= 2.0;
            Ey *= 2.0;
            Ez *= 2.0;
        } else {
            Ex += Db->vector_dot(intmats[0]);
            Ey += Db->vector_dot(intmats[1]);
            Ez += Db->vector_dot(intmats[2]);
        }
        field[3*n]   = Ex;
        field[3*n+1] = Ey;
        field[3*n+2] = Ez;
    }
    return EFP_RESULT_SUCCESS;
}

/*
 * Resetting of EFP options
 */
void EFP::set_options()
{
    enum efp_result res;

    molecule_ = Process::environment.molecule();

    struct efp_opts opts;
    memset(&opts, 0, sizeof(struct efp_opts));

    elst_enabled_ = options_.get_bool("EFP_ELST");
    pol_enabled_ = options_.get_bool("EFP_POL");
    disp_enabled_ = options_.get_bool("EFP_DISP");
    exch_enabled_ = options_.get_bool("EFP_EXCH");

    if (elst_enabled_)
        opts.terms |= EFP_TERM_ELEC;
    if (pol_enabled_)
        opts.terms |= EFP_TERM_POL;
    if (disp_enabled_)
        opts.terms |= EFP_TERM_DISP;
    if (exch_enabled_)
        opts.terms |= EFP_TERM_XR;
      //opts.terms |= EFP_TERM_CHTR;    // may be enabled in a future libefp release

    do_qm_ = options_.get_bool("QMEFP");
    qm_elst_enabled_ = do_qm_ && options_.get_bool("QMEFP_ELST");
    qm_pol_enabled_ = do_qm_ && options_.get_bool("QMEFP_POL");

    if (qm_elst_enabled_)
        opts.terms |= EFP_TERM_AI_ELEC;
    if (qm_pol_enabled_)
        opts.terms |= EFP_TERM_AI_POL;
      //opts.terms |= EFP_TERM_AI_DISP;  // may be enabled in future libefp release
      //opts.terms |= EFP_TERM_AI_XR;    // may be enabled in future libefp release
      //opts.terms |= EFP_TERM_AI_CHTR;  // may be enabled in future libefp release

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

    if ((res = efp_set_opts(efp_, &opts)))
        throw PsiException("EFP::set_options():efp_set_opts() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    do_grad_ = (options_.get_str("DERTYPE") == "FIRST");

    outfile->Printf("\n\n");
    outfile->Printf("%s", efp_banner());
    outfile->Printf("\n\n");

//**//    //outfile->Printf("  ==> Calculation Information (This section going away) <==\n\n");
//**//    //outfile->Printf("  Electrostatics damping: %12s\n", elst_damping_.c_str());
//**//
//**//    //if (disp_enabled_)
//**//    //    outfile->Printf("  Dispersion damping:     %12s\n", disp_damping_.c_str());
//**//
//**//    //outfile->Printf("  Electrostatics enabled:  %12d\n", elst_enabled_);
//**//    //outfile->Printf("  Polarization enabled:    %12d\n", pol_enabled_);
//**//    //outfile->Printf("  Dispersion enabled:      %12d\n", disp_enabled_);
//**//    //outfile->Printf("  Exchange enabled:        %12d\n", exch_enabled_);
//**//    //outfile->Printf("  Gradient enabled:        %12d\n", do_grad_);
//**//    //outfile->Printf("\n");
//**//
//**//    //outfile->Printf("\n");

    // sets call-back function to provide electric field from electrons
    if ((efp_set_electron_density_field_fn(efp_, electron_density_field_fn)))
        throw PsiException("EFP::set_options():efp_set_electron_density_field_fn() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

}

/*
 * Provide list of coordinates of quantum mechanical atoms to efp_set_point_charges
 */
void EFP::set_qm_atoms() {
    enum efp_result res;

    int natom = molecule_->natom();
    boost::shared_ptr<Vector> q (new Vector(natom));
    double * q_p = q->pointer();
    boost::shared_ptr<Vector> xyz (new Vector(3*natom));
    double * xyz_p = xyz->pointer();

    for (int A=0; A<natom; A++) {
        if (molecule_->Z(A) == 0.0)
            continue;
        q_p[A]       = molecule_->Z(A);
        xyz_p[3*A]   = molecule_->x(A);
        xyz_p[3*A+1] = molecule_->y(A);
        xyz_p[3*A+2] = molecule_->z(A);
    }

    if ((res = efp_set_point_charges(efp_, natom, q_p, xyz_p)))
        throw PsiException("EFP::set_qm_atoms():efp_set_point_charges() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
}

/*
 * Returns shared matrix containing the EFP contribution to the potential
 * felt by QM atoms, due to permanent EFP moments, in a SCF procedure.
 */
boost::shared_ptr<Matrix> EFP::modify_Fock_permanent()
{
    enum efp_result res;

    // get number of multipoles (charges, dipoles, quadrupoles, octupoles)
    size_t n_multipole = 0;
    if ((res = efp_get_multipole_count(efp_,&n_multipole)))
        throw PsiException("EFP::modify_Fock_permanent():efp_get_multipole_count(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    // multipole coordinates are stored array xyz.
    boost::shared_ptr<Vector> xyz  (new Vector(3*n_multipole));
    boost::shared_ptr<Vector> mult (new Vector((1+3+6+10)*n_multipole));

    // get multipoles from libefp
    //     dipoles stored as     x,y,z
    //     quadrupoles stored as xx,yy,zz,xy,xz,yz
    //     octupoles stored as   xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
    if ((res = efp_get_multipole_coordinates(efp_,xyz->pointer())))
        throw PsiException("EFP::modify_Fock_permanent():efp_get_multipole_coordinates(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
    if ((res = efp_get_multipole_values(efp_,mult->pointer())))
        throw PsiException("EFP::modify_Fock_permanent():efp_get_multipole_values(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    // Scale multipole integrals by multipole magnitudes.  The result goes into V
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<OneBodyAOInt> efp_ints(wfn->integral()->ao_efp_multipole_potential());

                               // 0    X    Y    Z      XX       YY       ZZ       XY       XZ       YZ
    const double prefacs[20] = { 1.0, 1.0, 1.0, 1.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0,
    //   XXX       YYY       ZZZ       XXY       XXZ       XYY       YYZ       XZZ       YZZ       XYZ
      1.0/15.0, 1.0/15.0, 1.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 3.0/15.0, 6.0/15.0};

    int nao = wfn->basisset()->nao();
    std::vector<SharedMatrix> mats;
    for(int i=0; i<20; ++i) {
        mats.push_back(SharedMatrix(new Matrix(nao, nao)));
    }

    // Cartesian basis one-electron EFP perturbation
    SharedMatrix V2(new Matrix("EFP permanent moment contribution to the Fock Matrix", nao, nao));

    // multipole contributions to Fock matrix
    double * xyz_p  = xyz->pointer();
    double * mult_p = mult->pointer();

    size_t max_natom = 0;
    for (int frag = 0; frag < nfrag_; frag++) {
        size_t natom = 0;
        if ((res = efp_get_frag_atom_count(efp_,frag,&natom)))
            throw PsiException("EFP::modify_Fock_permanent():efp_get_frag_atom_count(): " +
                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
        if ( natom > max_natom ) max_natom = natom;
    }
    efp_atom * atoms = (efp_atom*)malloc(max_natom*sizeof(efp_atom));

    for (size_t n=0; n<n_multipole; n++) {
        for (int i=0; i<20; ++i) {
           mats[i]->zero();
        }
        Vector3 coords(xyz_p[n*3],xyz_p[n*3+1],xyz_p[n*3+2]);
        efp_ints->set_origin(coords);
        efp_ints->compute(mats);

        // add point charges from atoms to multipoles at atom center
        for (int frag=0; frag<nfrag_; frag++) {
            size_t natom = 0;
            if ((res = efp_get_frag_atom_count(efp_, frag, &natom)))
                throw PsiException("EFP::modify_Fock_permanent():efp_get_frag_atom_count(): " +
                    std::string (efp_result_to_string(res)),__FILE__,__LINE__);
            if ((res = efp_get_frag_atoms(efp_, frag, natom, atoms)))
                throw PsiException("EFP::modify_Fock_permanent():efp_get_frag_atoms(): " +
                    std::string (efp_result_to_string(res)),__FILE__,__LINE__);

            for (size_t i=0; i<natom; i++) {
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

        for (int i=0; i<20; ++i) {
            mats[i]->scale( -prefacs[i] * mult_p[20*n+i] );
            V2->add(mats[i]);
        }
    }
    free(atoms);

    boost::shared_ptr<PetiteList> pet(new PetiteList(wfn->basisset(),wfn->integral(),true));
    boost::shared_ptr<Matrix> U = pet->aotoso();

    boost::shared_ptr<Matrix> V = Matrix::triplet(U,V2,U,true,false,false);
    V->set_name("EFP permanent moment contribution to the Fock Matrix");

    return V;
}

/*
 * Returns shared matrix containing the EFP contribution to the potential
 * felt by QM atoms, due to EFP induced dipoles, in a SCF procedure.
 */
boost::shared_ptr<Matrix> EFP::modify_Fock_induced()
{
    enum efp_result res;

    // induced dipoles
    size_t n_id = 0;
    if ((res = efp_get_induced_dipole_count(efp_, &n_id)))
        throw PsiException("EFP::modify_Fock_induced():efp_get_induced_dipole_count() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    boost::shared_ptr<Vector> xyz_id (new Vector(3*n_id));
    if ((res = efp_get_induced_dipole_coordinates(efp_, xyz_id->pointer())))
        throw PsiException("EFP::modify_Fock_induced():efp_get_induced_dipole_coordinates() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    boost::shared_ptr<Vector> id (new Vector(3*n_id));
    if ((res = efp_get_induced_dipole_values(efp_, id->pointer())))
        throw PsiException("EFP::modify_Fock_induced():efp_get_induced_dipole_values() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    boost::shared_ptr<Vector> idt (new Vector(3*n_id));
    if ((res = efp_get_induced_dipole_conj_values(efp_, idt->pointer())))
        throw PsiException("EFP::modify_Fock_induced():efp_get_induced_dipole_conj_values() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    // take average of induced dipole and conjugate
    id->add(idt);
    id->scale(0.5);

    // scale field integrals by induced dipole magnitudes.  the result goes into V
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<OneBodyAOInt> field_ints(wfn->integral()->electric_field());

    int nao = wfn->basisset()->nao();
    std::vector<SharedMatrix> mats;
    for (int i=0; i<3; ++i)
        mats.push_back(SharedMatrix(new Matrix(nao, nao)));

    // Cartesian basis one-electron EFP perturbation
    SharedMatrix V2(new Matrix("EFP induced dipole contribution to the Fock Matrix", nao, nao));

    // induced dipole contributions to Fock matrix
    // multipole contributions to Fock matrix
    double *xyz_p  = xyz_id->pointer();
    double *mult_p = id->pointer();
    for (size_t n=0; n<n_id; n++) {
        for (int i=0; i<3; ++i) {
           mats[i]->zero();
        }
        Vector3 coords(xyz_p[n*3],xyz_p[n*3+1],xyz_p[n*3+2]);
        field_ints->set_origin(coords);
        field_ints->compute(mats);
        // only dealing with dipoles here:
        for (int i=0; i<3; ++i) {
            mats[i]->scale(-mult_p[3*n+i]);
            V2->add(mats[i]);
        }
    }

    boost::shared_ptr<PetiteList> pet(new PetiteList(wfn->basisset(),wfn->integral(),true));
    boost::shared_ptr<Matrix> U = pet->aotoso();
    boost::shared_ptr<Matrix> V = Matrix::triplet(U,V2,U,true,false,false);
    V->set_name("EFP induced dipole contribution to the Fock Matrix");

    return V;
}

/*
 * compute efp contribution to scf energy
 */
double EFP::scf_energy_update()
{
    enum efp_result res;
    double efp_energy;

    if ((res = efp_get_wavefunction_dependent_energy(efp_, &efp_energy)))
        throw PsiException("EFP::scf_energy_update(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    return efp_energy;
}

/*
 * Compute efp energy components and/or gradient
 */
void EFP::compute() {
    enum efp_result res;
    struct efp_energy energy;

    if ((res = efp_compute(efp_, do_grad_ ? 1 : 0)))
        throw PsiException("EFP::compute():efp_compute() " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    if ((res = efp_get_energy(efp_, &energy)))
        throw PsiException("EFP::compute():efp_get_energy(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    if (do_grad_) {
        //SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag_, 6));
        //double ** psmgrad = smgrad->pointer();
        //if ((res = efp_get_gradient(efp_, psmgrad[0])))
        //    throw PsiException("EFP::compute():efp_get_gradient(): " +
        //        std::string (efp_result_to_string(res)),__FILE__,__LINE__);

        //outfile->Printf("  ==> EFP Gradient <==\n\n");

        //for (int i=0; i<nfrag_; i++) {
        //    for (int j=0; j<6; j++) {
        //        outfile->Printf("%14.6lf", psmgrad[i][j]);
        //    }
        //    outfile->Printf("\n");
        //}
        //outfile->Printf("\n");

        //torque_ = smgrad;
        //psi::Process::environment.set_efp_torque(smgrad);
        //smgrad->print();
    }

    outfile->Printf("  ==> Energetics <==\n\n");

    outfile->Printf("  EFP/EFP Electrostatics Energy = %20.12f [H] %s\n",
        energy.electrostatic + energy.charge_penetration, elst_enabled_ ? "*" : "");
    if (do_qm_) {
        outfile->Printf("  QM-Nuc/EFP Electrostatics Energy = %20.12f [H] %s\n",
            energy.electrostatic_point_charges, qm_elst_enabled_ ? "*" : "");
        outfile->Printf("  Total Electrostatics Energy = %20.12f [H] %s\n",
            energy.electrostatic + energy.charge_penetration + energy.electrostatic_point_charges,
            (elst_enabled_ || qm_elst_enabled_) ? "*" : "");
    }
    outfile->Printf("  %7s Polarization Energy =   %20.12f [H] %s\n", do_qm_ ? "" : "EFP/EFP",
        energy.polarization, pol_enabled_ ? "*" : "");
    outfile->Printf("  EFP/EFP Dispersion Energy =     %20.12f [H] %s\n",
        energy.dispersion, disp_enabled_ ? "*" : "");
    outfile->Printf("  EFP/EFP Exchange Energy =       %20.12f [H] %s\n",
        energy.exchange_repulsion, exch_enabled_ ? "*" : "");
    outfile->Printf("  Total Energy =                  %20.12f [H] %s\n",
        energy.total, "*");

    Process::environment.globals["EFP ELST ENERGY"] = energy.electrostatic +
                                                      energy.charge_penetration +
                                                      energy.electrostatic_point_charges;
    Process::environment.globals["EFP IND ENERGY"] = energy.polarization;
    Process::environment.globals["EFP DISP ENERGY"] = energy.dispersion;
    Process::environment.globals["EFP EXCH ENERGY"] = energy.exchange_repulsion;
    Process::environment.globals["EFP TOTAL ENERGY"] = energy.total;
    Process::environment.globals["CURRENT ENERGY"] = energy.total;
}

/*
 * Get number of fragments
 */
int EFP::get_frag_count(void)
{
    enum efp_result res;
    size_t n=0;

    if ((res = efp_get_frag_count(efp_, &n)))
        throw PsiException("EFP::get_frag_count(): " +
            std::string (efp_result_to_string(res)),__FILE__,__LINE__);

    nfrag_ = n;
    return n;
}

/*
 * Prints efp geometry is styly of Molecule
 */
void EFP::print_efp_geometry()
{
    if ( nfrag_ == 0 ) return;

    Molecule::GeometryUnits units = Process::environment.molecule()->units();
    bool is_angstrom = ( units == Molecule::Angstrom ) ;

    outfile->Printf("\n");
    outfile->Printf("  ==> EFP Geometry <==\n\n");
    outfile->Printf("    Geometry (in %s):\n\n", //, charge = %d, multiplicity = %d:\n\n",
        is_angstrom ? "Angstrom" : "Bohr"); //, get_frag_multiplicity, multiplicity_);
    outfile->Printf("       Center              X                  Y                   Z       \n");
    outfile->Printf("    ------------   -----------------  -----------------  -----------------\n");

    for (int frag=0; frag<nfrag_; frag++) {
        double * xyz = get_frag_atom_coord(frag);
        std::vector<std::string> symbol = get_frag_atom_label(frag);
        int natom = symbol.size();
        for (int i=0; i<natom; i++) {
            regex_match(symbol[i], reMatches, efpAtomSymbol);
            outfile->Printf( "    %8s%4s ", reMatches[1].str().c_str(), "");
            for (int j=0; j<3; j++) {
                outfile->Printf("  %17.12lf", xyz[i*3+j] * ( is_angstrom ? pc_bohr2angstroms : 1.0 ) );
            }
            outfile->Printf(" (EFP %3i)\n", frag+1);
        }
    }
    outfile->Printf("\n");
}

/*
 * Print simple private members of class
 */
void EFP::print_out() {

    outfile->Printf("  ==> EFP/EFP Setup <==\n\n");
    outfile->Printf("  Number of EFP fragments: %12d\n", nfrag_);
    outfile->Printf("  Electrostatics enabled?: %12s\n", elst_enabled_ ? "true" : "false");
    outfile->Printf("  Polarization enabled?:   %12s\n", pol_enabled_ ? "true" : "false");
    outfile->Printf("  Dispersion enabled?:     %12s\n", disp_enabled_ ? "true" : "false");
    outfile->Printf("  Exchange enabled?:       %12s\n", exch_enabled_ ? "true" : "false");
    if (elst_enabled_)
        outfile->Printf("  Electrostatics damping:  %12s\n", elst_damping_.c_str());
    if (pol_enabled_)
        outfile->Printf("  Polarization damping:    %12s\n", pol_damping_.c_str());
    if (disp_enabled_)
        outfile->Printf("  Dispersion damping:      %12s\n", disp_damping_.c_str());
    outfile->Printf("  Gradient enabled?:       %12s\n", do_grad_ ? "true" : "false");

    if (do_qm_) {
        outfile->Printf("\n  ==> QM/EFP Setup <==\n\n");
        outfile->Printf("  Number of QM fragments:  %12d\n", molecule_->nfragments());
        outfile->Printf("  Electrostatics enabled?: %12s\n", qm_elst_enabled_ ? "true" : "false");
        outfile->Printf("  Polarization enabled?:   %12s\n", qm_pol_enabled_ ? "true" : "false");
        outfile->Printf("  Dispersion enabled?:     %12s\n", "undefined");
        outfile->Printf("  Exchange enabled?:       %12s\n", "undefined");
    }

    print_efp_geometry();

    if (do_qm_) {
        outfile->Printf("  ==> QM Geometry <==\n\n");
        molecule_->print();
    }
}

} } // End namespaces


///*
// * Get COM for a fragment from libefp coordinates
// */
//double * EFP::get_com(int frag_idx)
//{
//    if (frag_idx >= nfrag_)
//        return NULL;
//
//    double * xyzabc = new double [6*nfrag_];
//    efp_get_coordinates(efp_, xyzabc);
//
//    double * com = new double[3];
//    com[0] = xyzabc[6*frag_idx+0];
//    com[1] = xyzabc[6*frag_idx+1];
//    com[2] = xyzabc[6*frag_idx+2];
//
//    return com;
//}

///*
// * Returns number of EFP atoms
// */
//int EFP::efp_natom()
//{
//    if (nfrag_ == 0)
//        return 0;
//    int natom = 0;
//    for (int frag=0; frag<nfrag_; frag++) {
//        std::vector<std::string> symbol = get_frag_atom_label(frag);
//        natom += symbol.size();
//    }
//    return natom;
//}

//boost::shared_ptr<Matrix> EFP::EFP_nuclear_potential()
//{
//    enum efp_result res;
//
//    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
//
//    int nbf = wfn->basisset()->nbf();
//
//    boost::shared_ptr<Matrix> V (new Matrix(nbf,nbf));
//
//    for (int frag=0; frag<nfrag_; frag++) {
//        size_t natom = 0;
//        if ((res = efp_get_frag_atom_count(efp_, frag, &natom)))
//            throw PsiException("EFP::EFP_nuclear_potential():efp_get_frag_atom_count() " +
//                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//        efp_atom * atoms = (efp_atom*)malloc(natom*sizeof(efp_atom));
//        if ((res = efp_get_frag_atoms(efp_, frag, natom, atoms)))
//            throw PsiException("EFP::EFP_nuclear_potential():efp_get_frag_atoms() " +
//                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//        SharedMatrix V_charge(new Matrix("External Potential (Charges)", nbf, nbf));
//
//        SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", natom, 4));
//        double** Zxyzp = Zxyz->pointer();
//        for (int i=0; i<natom; i++) {
//            Zxyzp[i][0] = atoms[i].znuc;
//            Zxyzp[i][1] = atoms[i].x;
//            Zxyzp[i][2] = atoms[i].y;
//            Zxyzp[i][3] = atoms[i].z;
//        }
//
//        boost::shared_ptr<PotentialInt> pot(static_cast<PotentialInt*>(wfn->integral()->ao_potential()));
//        pot->set_charge_field(Zxyz);
//        pot->compute(V_charge);
//
//        // add to efp fock contribution
//        V->add(V_charge);
//
//        V_charge.reset();
//        pot.reset();
//    }
//    return V;
//}

///*
// * Get gradient of the interaction energy of the EFP electrostatics with the QM nuclei (point charges)
// */
//boost::shared_ptr<Vector> EFP::get_electrostatic_gradient()
//{
//    enum efp_result res;
//
//    int natom = molecule_->natom();
//    boost::shared_ptr<Vector> grad ( new Vector( 3*natom ) );
//    double * grad_p = grad->pointer();
//
//    // verify natom matches the number of point charges in efp
//    size_t n_ptc;
//    if ((res = efp_get_point_charge_count(efp_, &n_ptc)))
//        throw PsiException("EFP::get_electrostatic_gradient():efp_get_point_charge_count() " +
//            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//    if (n_ptc != natom)
//        throw PsiException("EFP::get_electrostatic_gradient(): natom does not match number of point charges in efp",
//            __FILE__,__LINE__);
//
//    if ((res = efp_get_point_charge_gradient(efp_, grad_p)))
//        throw PsiException("EFP::get_electrostatic_gradient():efp_get_point_charge_gradient() " +
//            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//    return grad;
//}

///*
// * Add fragment by name
// */
//void EFP::add_fragment(std::string fname)
//{
//    enum efp_result res;
//
//    if ((res = efp_add_fragment(efp_, fname.c_str())))
//        throw PsiException("EFP::add_fragment(): " + fname +
//            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//}

///*
// * Set points or xyzabc coordinates for all fragments simultaneously
// */
//void EFP::set_coordinates(int type, double * coords)
//{
//    enum efp_result res;
//    enum efp_coord_type ctype;
//
//    if(type == 0)
//        ctype = EFP_COORD_TYPE_XYZABC;
//    else if(type == 1)
//        ctype = EFP_COORD_TYPE_POINTS;
//    else if(type == 2)
//        ctype = EFP_COORD_TYPE_ROTMAT;
//
//    if ((res = efp_set_coordinates(efp_, ctype, coords)))
//        throw PsiException("EFP::set_coordinates(): " +
//            std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//}

///*
// * Compute QM-Nuc/EFP-Nuc NRE
// */
//double EFP::EFP_QM_nuclear_repulsion_energy()
//{
//    double nu = 0.0;
//    boost::shared_ptr<Molecule> mol = Process::environment.molecule();
//
//    for (int frag=0; frag<nfrag_; frag++) {
//        size_t natom = 0;
//        if ((res = efp_get_frag_atom_count(efp_, frag, &natom)))
//            throw PsiException("EFP::EFP_QM_nuclear_repulsion_energy():efp_frag_atom_count() " +
//                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//        efp_atom * atoms = (efp_atom*)malloc(natom*sizeof(efp_atom));
//        if ((res = efp_get_frag_atoms(efp_, frag, natom, atoms)))
//            throw PsiException("EFP::EFP_QM_nuclear_repulsion_energy():efp_get_frag_atoms() " +
//                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
//
//        for (int i=0; i<natom; i++) {
//            double znuc = atoms[i].znuc;
//            double x    = atoms[i].x;
//            double y    = atoms[i].y;
//            double z    = atoms[i].z;
//
//            for (int j=0; j<mol->natom(); j++) {
//                double dx = x - mol->x(j);
//                double dy = y - mol->y(j);
//                double dz = z - mol->z(j);
//                double r  = sqrt(dx*dx+dy*dy+dz*dz);
//                nu += znuc * mol->Z(j) / r;
//            }
//        }
//    }
//    return nu;
//}

