/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*
 *  sad.cc
 *
 * Routines for the high-maintenance SAD guess
 * and dual-basis projections
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>
#include <libfock/jk.h>

#include "hf.h"
#include "sad.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

SADGuess::SADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options) :
    basis_(basis), nalpha_(nalpha), nbeta_(nbeta), options_(options)
{
    common_init();
}
SADGuess::~SADGuess()
{
}
void SADGuess::common_init()
{
    molecule_ = basis_->molecule();

    boost::shared_ptr<IntegralFactory> ints(new IntegralFactory(basis_));
    boost::shared_ptr<PetiteList> petite(new PetiteList(basis_,ints));
    AO2SO_ =  petite->aotoso();

    print_ = options_.get_int("SAD_PRINT");
    debug_ = options_.get_int("DEBUG");
    if(options_["SOCC"].size()>0||options_["DOCC"].size()>0)
       PSIEXCEPTION("SAD guess not implemented for user-specified SOCCs and/or DOCCs yet");
}
void SADGuess::compute_guess()
{
    
    timer_on("SAD Guess");
    form_D();
    form_C();
    timer_off("SAD Guess");
}
void SADGuess::form_D()
{
    // Build Neutral D in AO basis (block diagonal)
    SharedMatrix DAO = form_D_AO();

    // Transform Neutral D from AO to SO basis
    Da_ = SharedMatrix(new Matrix("Da SAD",AO2SO_->colspi(),AO2SO_->colspi()));

    double* temp = new double[AO2SO_->rowspi()[0] * (ULI) AO2SO_->max_ncol()];
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nao = AO2SO_->rowspi()[h];
        int nso = AO2SO_->colspi()[h];
        if (!nao || !nso) continue;

        double** DAOp = DAO->pointer();
        double** DSOp = Da_->pointer(h);
        double** Up = AO2SO_->pointer(h);

        C_DGEMM('N','N',nao,nso,nao,1.0,DAOp[0],nao,Up[0],nso,0.0,temp,nso);
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,temp,nso,0.0,DSOp[0],nso);
    }
    delete[] temp;

    // Scale Da to true electron count
    double npair = 0.0;
    for (int A = 0; A < molecule_->natom(); A++) {
        npair += 0.5 * molecule_->Z(A);
    }
    Da_->scale(((double) nalpha_) / npair);

    // Build/Scale Db if needed
    if (nalpha_ == nbeta_) {
        Db_ = Da_;
    } else {
        Db_ = SharedMatrix(Da_->clone());
        Db_->set_name("Db SAD");
        Db_->scale(((double) nbeta_) / ((double) nalpha_));
    }

    if (debug_) {
        Da_->print();
        Db_->print();
    }
}
void SADGuess::form_C()
{
    Ca_ = Da_->partial_cholesky_factorize(options_.get_double("SAD_CHOL_TOLERANCE"));
    Ca_->set_name("Ca SAD");
    if (nalpha_ == nbeta_) {

        Cb_ = Ca_;
    } else {
        Cb_ = SharedMatrix(Ca_->clone());
        Cb_->set_name("Cb SAD");
        Cb_->scale(sqrt(((double)nbeta_)/((double)nalpha_)));
    }

    if (debug_) {
        Ca_->print();
        Cb_->print();
    }
}
SharedMatrix SADGuess::form_D_AO()
{
    std::vector<boost::shared_ptr<BasisSet> > atomic_bases;

    if (print_ > 6) {
        outfile->Printf("\n  Constructing atomic basis sets\n  Molecule:\n");
        molecule_->print();
    }

    //Build the atomic basis sets for libmints use in UHF
    for (int A = 0; A<molecule_->natom(); A++) {
        std::stringstream mol_string;
        mol_string << std::endl << basis_->molecule()->label(A) << std::endl;
        boost::shared_ptr<Molecule> atom_mol = Molecule::create_molecule_from_string(mol_string.str());
        atom_mol->reset_point_group("C1"); // Booo symmetry

        boost::shared_ptr<BasisSet> atom_bas = BasisSet::pyconstruct_orbital(atom_mol, "BASIS",
                                                basis_->molecule()->atom_entry(A)->basisset());
        atomic_bases.push_back(atom_bas);

        if (print_ > 6) {
            outfile->Printf("  SAD: Atomic Basis Set %d\n", A);
            atomic_bases[A]->molecule()->print();
            outfile->Printf("\n");
            atomic_bases[A]->print("outfile");
            outfile->Printf("\n");
        }
    }

    //Spin occupations per atom, to be determined by Hund's Rules
    //or user input
    std::vector<int> nalpha(molecule_->natom(), 0);
    std::vector<int> nbeta(molecule_->natom(), 0);
    std::vector<int> nelec(molecule_->natom(), 0);
    std::vector<int> nhigh(molecule_->natom(), 0);
    int tot_elec = 0;

    //Ground state high spin occupency array, atoms 0 to 36 (see Giffith's Quantum Mechanics, pp. 217)
    //For 37 to 86, save for f-block: Atomic, Molecular, & Optical Physics Handbook, Ed. Gordon W. F. Drake, American Institute of Physics, Woodbury, New York, USA, 1996.
    const int reference_S[] = {0,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,6,5,4,3,2,1,0,1,2,3,2,1,0,1,0,1,2,5,6,5,4,3,0,1,0,1,2,3,2,1,0,1,0,1,0,3,4,5,6,7,8,5,4,3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,2,1,0};
    const int MAX_Z = 86;

    if (print_ > 1)
        outfile->Printf("  Determining Atomic Occupations\n");

    for (int A = 0; A<molecule_->natom(); A++) {
        int Z = molecule_->Z(A);
        if (Z>MAX_Z) {
            throw std::domain_error(" Only Atoms up to 86 (Rn) are currently supported with SAD Guess");
        }
        nhigh[A] = reference_S[Z];
        nelec[A] = Z;
        tot_elec+= nelec[A];
        nbeta[A] = (nelec[A]-nhigh[A])/2;
        nalpha[A] = nelec[A]-nbeta[A];
        if (print_ > 1)
            outfile->Printf("  Atom %d, Z = %d, nelec = %d, nhigh = %d, nalpha = %d, nbeta = %d\n",A,Z,nelec[A],nhigh[A],nalpha[A],nbeta[A]);
    }


    // Determine redundant atoms
    std::vector<int> unique_indices(molecule_->natom(), 0); // All atoms to representative unique atom
    std::vector<int> atomic_indices(molecule_->natom(), 0); // unique atom to first representative atom
    std::vector<int> offset_indices(molecule_->natom(), 0); // unique atom index to rank
    int nunique = 0;
    for (int l = 0; l < molecule_->natom(); l++) {
        unique_indices[l] = l;
        atomic_indices[l] = l;
    }

    // This is all overkill for now, but lets leave it incase we do something oddball in the future
    for (int l = 0; l < molecule_->natom() - 1; l++) {
        for (int m = l + 1; m < molecule_->natom(); m++) {
            if (unique_indices[m] != m)
                continue; //Already assigned
            if (molecule_->Z(l) != molecule_->Z(m))
                continue;
            if (nalpha[l] != nalpha[m])
                continue;
            if (nbeta[l] != nbeta[m])
                continue;
            if (nhigh[l] != nhigh[m])
                continue;
            if (nelec[l] != nelec[m])
                continue;
            if (atomic_bases[l]->nbf() != atomic_bases[m]->nbf())
                continue;
            if (atomic_bases[l]->nshell() != atomic_bases[m]->nshell())
                continue;
            if (atomic_bases[l]->nprimitive() != atomic_bases[m]->nprimitive())
                continue;
            if (atomic_bases[l]->max_am() != atomic_bases[m]->max_am())
                continue;
            if (atomic_bases[l]->max_nprimitive() != atomic_bases[m]->max_nprimitive())
                continue;
            if (atomic_bases[l]->has_puream() !=  atomic_bases[m]->has_puream())
                continue;

            // Semi-Rigorous match obtained
            unique_indices[m] = l;
        }
    }
    for (int l = 0; l < molecule_->natom(); l++) {
        if (unique_indices[l] == l) {
            atomic_indices[nunique] = l;
            offset_indices[l] = nunique;
            nunique++;
        }
    }

    //Atomic D matrices within the atom specific AO basis
    std::vector<SharedMatrix> atomic_D;
    for (int A = 0; A<nunique; A++) {
        int nbf = atomic_bases[atomic_indices[A]]->nbf();
        SharedMatrix dtmp(new Matrix("Atomic D", nbf, nbf));
        atomic_D.push_back(dtmp);
    }

    if (print_ > 1)
        outfile->Printf("\n  Performing Atomic UHF Computations:\n");
    for (int A = 0; A<nunique; A++) {
        int index = atomic_indices[A];
        if (print_ > 1)
            outfile->Printf("\n  UHF Computation for Unique Atom %d which is Atom %d:",A, index);
        get_uhf_atomic_density(atomic_bases[index], nelec[index], nhigh[index], atomic_D[A]);
        if (print_ > 1)
            outfile->Printf("Finished UHF Computation!\n");
    }
    if (print_)
        outfile->Printf("\n");



    //Add atomic_D into D (scale by 1/2, we like effective pairs)
    SharedMatrix DAO = SharedMatrix(new Matrix("D_SAD (AO)", basis_->nbf(), basis_->nbf()));
    for (int A = 0, offset = 0; A < molecule_->natom(); A++) {
        int norbs = atomic_bases[A]->nbf();
        int back_index = unique_indices[A];
        for (int m = 0; m<norbs; m++)
            for (int n = 0; n<norbs; n++)
                DAO->set(0, m+offset, n+offset,
                         0.5 * atomic_D[offset_indices[back_index]]->get(m,  n));
        offset += norbs;
    }

    if (debug_) {
        DAO->print();
    }

    return DAO;
}
void SADGuess::get_uhf_atomic_density(boost::shared_ptr<BasisSet> bas, int nelec, int nhigh, SharedMatrix D)
{
    boost::shared_ptr<Molecule> mol = bas->molecule();
    mol->update_geometry();
    if (print_ > 1){
        mol->print();
    }

    int nbeta = (nelec - nhigh) / 2;
    int nalpha = nelec - nbeta;
    int natom = mol->natom();
    int norbs = bas->nbf();
    int Z = bas->molecule()->Z(0);

    if (nalpha > norbs || nbeta > norbs) throw PSIEXCEPTION("Atom has more electrons than basis functions.");

    if (print_ > 1) {
        outfile->Printf("\n");
        bas->print("outfile");
        outfile->Printf("  Occupation: nalpha = %d, nbeta = %d, norbs = %d\n",nalpha,nbeta,norbs);
        outfile->Printf("\n  Atom:\n");
        mol->print();
    }

    if (natom != 1) {
        throw std::domain_error("SAD Atomic UHF has been given a molecule, not an atom");
    }

    IntegralFactory integral(bas, bas, bas, bas);
    MatrixFactory mat;
    mat.init_with(1,&norbs,&norbs);
    OneBodyAOInt *S_ints = integral.ao_overlap();
    OneBodyAOInt *T_ints = integral.ao_kinetic();
    OneBodyAOInt *V_ints = integral.ao_potential();

    // Compute overlap S and orthogonalizer X;
    SharedMatrix S(mat.create_matrix("Overlap Matrix"));
    S_ints->compute(S);

    SharedMatrix X = S->clone();
    X->power(-0.5, 1.e-10);
    X->set_name("Orthogonalizer X^-1/2 Matrix");

    if (print_ > 6) {
        S->print();
        X->print();
    }

    //Compute H
    SharedMatrix T(mat.create_matrix("T"));
    T_ints->compute(T);
    SharedMatrix V(mat.create_matrix("V"));
    V_ints->compute(V);
    SharedMatrix H(mat.create_matrix("Core Hamiltonian Matrix H"));
    H->zero();
    H->add(T);
    H->add(V);

    T.reset();
    V.reset();
    delete S_ints;
    delete T_ints;
    delete V_ints;

    if (print_ > 6) {
        H->print();
    }

    // Init temps
    SharedMatrix Ca(mat.create_matrix("Ca"));
    SharedMatrix Cb(mat.create_matrix("Cb"));

    SharedMatrix Da(mat.create_matrix("Da"));
    SharedMatrix Db(mat.create_matrix("Db"));

    SharedMatrix gradient_a(mat.create_matrix("gradient_a"));
    SharedMatrix gradient_b(mat.create_matrix("gradient_b"));

    SharedMatrix Fa(mat.create_matrix("Fa"));
    SharedMatrix Fb(mat.create_matrix("Fb"));

    // Factional occupation
    SharedVector occ_a, occ_b;
    if (options_.get_double("SAD_FRAC_OCC")){
        int nfzc = 0, nact = 0;
        if (Z <= 2){
            nfzc = 0;
            nact = 1;
        }
        else if (Z <= 4){
            nfzc = 1;
            nact = 1;
        }
        else if (Z <= 10){
            nfzc = 2;
            nact = 3;
        }
        else if (Z <= 18){
            nfzc = 5;
            nact = 4;
        }
        else if (Z <= 36){
            nfzc = 9;
            nact = 9;
        }
        else if (Z <= 54){
            nfzc = 18;
            nact = 9;
        }
        else if (Z <= 54){
            nfzc = 18;
            nact = 9;
        }
        else if (Z <= 86){
            nfzc = 27;
            nact = 16;
        }
        else{
            throw PSIEXCEPTION("SAD: Fractional occupations are not supported beyond Radeon");
        }

        nalpha = nfzc + nact;
        nbeta = nalpha;
        double frac_act = std::pow(((double)(Z - nfzc * 2)) / ((double)nact * 2), 0.5);

        occ_a = boost::shared_ptr<Vector>(new Vector("Alpha fractional occupation", nalpha));
        for (size_t x = 0; x<nfzc; x++) occ_a->set(x, 1.0);
        for (size_t x = nfzc; x<nalpha; x++) occ_a->set(x, frac_act);
        occ_b = boost::shared_ptr<Vector>(occ_a->clone());

    }
    else{

        // Conventional occupations
        occ_a = boost::shared_ptr<Vector>(new Vector("Alpha occupation", nalpha));
        for (size_t x = 0; x<nalpha; x++) occ_a->set(x, 1.0);
        occ_b = boost::shared_ptr<Vector>(new Vector("Beta occupation", nbeta));
        for (size_t x = 0; x<nbeta; x++) occ_b->set(x, 1.0);
    }

    SharedMatrix Ca_occ(new Matrix("Ca occupied", norbs, nalpha));
    SharedMatrix Cb_occ(new Matrix("Ca occupied", norbs, nbeta));

    //Compute initial Cx, Dx, and D from core guess
    form_C_and_D(nalpha, norbs, X, H, Ca, Ca_occ, occ_a, Da);
    form_C_and_D(nbeta, norbs, X, H, Cb, Cb_occ, occ_b, Db);

    D->zero();
    D->add(Da);
    D->add(Db);

    if (print_ > 6) {
        Ca->print();
        Cb->print();
        Ca_occ->print();
        Cb_occ->print();
        Da->print();
        Db->print();
        D->print();
    }

    //Compute inital E for reference
    double E = D->vector_dot(H);
    E *= 0.5;

    double E_tol = options_.get_double("SAD_E_CONVERGENCE");
    double D_tol = options_.get_double("SAD_D_CONVERGENCE");
    int maxiter = options_.get_int("SAD_MAXITER");

    double E_old = E;
    int iteration = 0;

    // Setup DIIS
    DIISManager diis_manager(6, "SAD DIIS", DIISManager::LargestError, DIISManager::InCore);
    diis_manager.set_error_vector_size(2, DIISEntry::Matrix, gradient_a.get(), DIISEntry::Matrix, gradient_b.get());
    diis_manager.set_vector_size(2, DIISEntry::Matrix, Fa.get(), DIISEntry::Matrix, Fb.get());

    // Setup JK
    std::unique_ptr<JK> jk;

    // Need a very special auxiliary basis here
    if (options_.get_str("SAD_SCF_TYPE") == "DF"){
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_orbital(bas->molecule(), "BASIS",
                                                options_.get_str("DF_BASIS_SAD"));

        DFJK* dfjk = new DFJK(bas, auxiliary);
        dfjk->set_unit(PSIF_SAD);
        if (options_["DF_INTS_NUM_THREADS"].has_changed())
            dfjk->set_df_ints_num_threads(options_.get_int("DF_INTS_NUM_THREADS"));
        jk = std::unique_ptr<JK>(dfjk);
    }
    else if (options_.get_str("SAD_SCF_TYPE") == "DIRECT"){
        DirectJK* directjk(new DirectJK(bas));
        if (options_["DF_INTS_NUM_THREADS"].has_changed())
            directjk->set_df_ints_num_threads(options_.get_int("DF_INTS_NUM_THREADS"));
        jk = std::unique_ptr<JK>(directjk);
    }
    else {
        std::stringstream msg;
        msg << "SAD: JK type of " << options_.get_str("SAD_SCF_TYPE") << " not understood.\n";
        throw PSIEXCEPTION(msg.str());
    }

    jk->set_memory((ULI)(0.5 * (Process::environment.get_memory() / 8L)));
    jk->initialize();
    if (print_ > 1)
        jk->print_header();

    // These are static so lets just grab them now
    std::vector<SharedMatrix> & jkC = jk->C_left();
    jkC.push_back(Ca_occ);
    jkC.push_back(Cb_occ);
    const std::vector<SharedMatrix> & Jvec = jk->J();
    const std::vector<SharedMatrix> & Kvec = jk->K();

    // Print a header
    bool converged = false;
    if (print_ > 1) {
        outfile->Printf( "\n  Initial Atomic UHF Energy:    %14.10f\n\n",E);
        outfile->Printf( "                                         Total Energy            Delta E              Density RMS\n\n");

    }

    // Run the iterations
    do {

        iteration++;

        // Copy the old values over for error analysis
        E_old = E;

        // Compute JK matrices
        jk->compute();

        // Form Fa and Fb
        Fa->copy(H);
        Fa->add(Jvec[0]);
        Fa->add(Jvec[1]);

        Fb->copy(Fa);

        Fa->subtract(Kvec[0]);
        Fb->subtract(Kvec[1]);

        // Compute E
        E  = H->vector_dot(D);
        E += Da->vector_dot(Fa);
        E += Db->vector_dot(Fb);
        E *= 0.5;

        double deltaE = fabs(E-E_old);

        // Build Gradient
        form_gradient(norbs, gradient_a, Fa, Da, S, X);
        form_gradient(norbs, gradient_b, Fb, Db, S, X);
        double Drms = 0.5 * (gradient_a->rms() + gradient_b->rms());

        // Add and extrapolate DIIS
        diis_manager.add_entry(4, gradient_a.get(), gradient_b.get(), Fa.get(), Fb.get());
        diis_manager.extrapolate(2, Fa.get(), Fb.get());

        //Diagonalize Fa and Fb to from Ca and Cb and Da and Db
        form_C_and_D(nalpha, norbs, X, Fa, Ca, Ca_occ, occ_a, Da);
        form_C_and_D(nbeta, norbs, X, Fb, Cb, Cb_occ, occ_b, Db);

        //Form D
        D->copy(Da);
        D->add(Db);

        if (print_ > 6) {
            H->print();
            Fa->print();
            Fb->print();
            Ca->print();
            Cb->print();
            Da->print();
            Db->print();
            D->print();
        }
        if (print_ > 1)
            outfile->Printf( "  @Atomic UHF iteration %3d energy: %20.14f    %20.14f %20.14f\n", iteration, E, E-E_old, Drms);

        //Check convergence
        if (iteration > 1 && deltaE < E_tol && Drms < D_tol)
            converged = true;

        if (iteration > maxiter) {
            outfile->Printf( "\n WARNING: Atomic UHF is not converging! Try casting from a smaller basis or call Rob at CCMST.\n");
            break;
        }

    } while (!converged);

    if (converged && print_ > 1)
        outfile->Printf( "  @Atomic UHF Final Energy for atom %s: %20.14f\n", mol->symbol(0).c_str(),E);

}
void SADGuess::form_gradient(int norbs, SharedMatrix grad, SharedMatrix F, SharedMatrix D,
                             SharedMatrix S, SharedMatrix X)
{
    SharedMatrix Scratch1(new Matrix("Scratch1", norbs, norbs));
    SharedMatrix Scratch2(new Matrix("Scratch2", norbs, norbs));

    // FDS
    Scratch1->gemm(false, false, 1.0, F, D, 0.0);
    Scratch2->gemm(false, false, 1.0, Scratch1, S, 0.0);

    // SDF
    Scratch1->copy(Scratch2);
    Scratch1->transpose_this();

    // FDS - SDF
    grad->copy(Scratch2);
    grad->subtract(Scratch1);

    // Level it out X(FDS - SDF)X
    Scratch1->gemm(false, false, 1.0, X, grad, 0.0);
    grad->gemm(false, false, 1.0, Scratch1, X, 0.0);

    Scratch1.reset();
    Scratch2.reset();
}


void SADGuess::form_C_and_D(int nocc, int norbs, SharedMatrix X, SharedMatrix F,
                                        SharedMatrix C, SharedMatrix Cocc, SharedVector occ,
                                        SharedMatrix D)
{
    if (nocc == 0) return;

    //Forms C in the AO basis for SAD Guesses
    SharedMatrix Scratch1(new Matrix("Scratch1", norbs, norbs));
    SharedMatrix Scratch2(new Matrix("Scratch2", norbs, norbs));

    // Form Fp = XFX
    Scratch1->gemm(true, false, 1.0, X, F, 0.0);
    Scratch2->gemm(false, false, 1.0, Scratch1, X, 0.0);

    SharedVector eigvals(new Vector("Eigenvalue scratch", norbs));
    Scratch2->diagonalize(Scratch1, eigvals);

    //Form C = XC'
    C->gemm(false, false, 1.0, X, Scratch1, 0.0);

    // Copy over Cocc
    double** Coccp = Cocc->pointer();
    double** Cp = C->pointer();
    for (int i = 0; i < norbs; i++){
        C_DCOPY(nocc, Cp[i], 1, Coccp[i], 1);
    }
    // Scale by occ
    for (int i = 0; i < nocc; i++){
        C_DSCAL(norbs, occ->get(i), &Cp[0][i], nocc);
    }
    //Form D = Cocc*Cocc'
    D->gemm(false, true, 1.0, Cocc, Cocc, 0.0);

    Scratch1.reset();
    Scratch2.reset();
}

void HF::compute_SAD_guess()
{

    boost::shared_ptr<SADGuess> guess(new SADGuess(basisset_, nalpha_, nbeta_, options_));
    guess->compute_guess();

    SharedMatrix Ca_sad = guess->Ca();
    SharedMatrix Cb_sad = guess->Cb();
    Da_->copy(guess->Da());
    Db_->copy(guess->Db());

    for (int h = 0; h < Da_->nirrep(); h++) {

        int nso = Ca_sad->rowspi()[h];
        int nmo = Ca_sad->colspi()[h];
        if (nmo > X_->colspi()[h])
            nmo = X_->colspi()[h];

        sad_nocc_[h] = nmo;

        if (!nso || !nmo) continue;

        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Ca2p = Ca_sad->pointer(h);
        double** Cb2p = Cb_sad->pointer(h);

        for (int i = 0; i < nso; i++) {
            ::memcpy((void*) Cap[i], (void*) Ca2p[i], nmo*sizeof(double));
            ::memcpy((void*) Cbp[i], (void*) Cb2p[i], nmo*sizeof(double));
        }
    }

    int temp_nocc;
    for (int h = 0 ; h < Da_->nirrep(); h++) {
        temp_nocc = sad_nocc_[h];
        sad_nocc_[h] = doccpi_[h];
        nalphapi_[h] = temp_nocc;
        nbetapi_[h]  = temp_nocc;
        doccpi_[h]   = temp_nocc;
        soccpi_[h]   = 0;
    }

    E_ = 0.0; // This is the -1th iteration
}
SharedMatrix HF::BasisProjection(SharedMatrix C_A, int* noccpi, boost::shared_ptr<BasisSet> old_basis, boost::shared_ptr<BasisSet> new_basis)
{

    //Based on Werner's method from Mol. Phys. 102, 21-22, 2311
    boost::shared_ptr<IntegralFactory> newfactory(new IntegralFactory(new_basis,new_basis,new_basis,new_basis));
    boost::shared_ptr<IntegralFactory> hybfactory(new IntegralFactory(old_basis,new_basis,old_basis,new_basis));
    boost::shared_ptr<OneBodySOInt> intBB(newfactory->so_overlap());
    boost::shared_ptr<OneBodySOInt> intAB(hybfactory->so_overlap());

    boost::shared_ptr<PetiteList> pet(new PetiteList(new_basis, newfactory));
    SharedMatrix AO2USO(pet->aotoso());

    SharedMatrix SAB(new Matrix("S_AB", C_A->nirrep(), C_A->rowspi(), AO2USO->colspi()));
    SharedMatrix SBB(new Matrix("S_BB", C_A->nirrep(), AO2USO->colspi(), AO2USO->colspi()));

    intAB->compute(SAB);
    intBB->compute(SBB);

    //SAB->print();
    //SBB->print();

    newfactory.reset();
    hybfactory.reset();
    intAB.reset();
    intBB.reset();
    pet.reset();

    // Constrained to the same symmetry at the moment, we can relax this soon
    SharedMatrix C_B(new Matrix("C_B", C_A->nirrep(), AO2USO->colspi(), noccpi));

    // Block over irreps (soon united irreps)
    for (int h = 0; h < C_A->nirrep(); h++) {

        int nocc = noccpi[h];
        int na = C_A->rowspi()[h];
        int nb = AO2USO->colspi()[h];

        if (nocc == 0 || na == 0 || nb == 0) continue;

        double** Ca = C_A->pointer(h);
        double** Cb = C_B->pointer(h);
        double** Sab = SAB->pointer(h);
        double** Sbb = SBB->pointer(h);

        int CholError = C_DPOTRF('L',nb,Sbb[0],nb);
        if (CholError !=0 )
            throw std::domain_error("S_BB Matrix Cholesky failed!");

        //Inversion (in place)
        int IError = C_DPOTRI('L',nb,Sbb[0],nb);
        if (IError !=0 )
            throw std::domain_error("S_BB Inversion Failed!");

        //LAPACK is smart and all, only uses half of the thing
        for (int m = 0; m<nb; m++)
            for (int n = 0; n<m; n++)
                Sbb[m][n] = Sbb[n][m];

        //Form T
        double** Temp1 = block_matrix(nb,nocc);
        C_DGEMM('T','N',nb,nocc,na,1.0,Sab[0],nb,Ca[0],nocc,0.0,Temp1[0],nocc);

        //outfile->Printf(" Temp1:\n");
        //print_mat(Temp1,nb,nocc,outfile);

        double** Temp2 = block_matrix(nb,nocc);
        C_DGEMM('N','N',nb,nocc,nb,1.0,Sbb[0],nb,Temp1[0],nocc,0.0,Temp2[0],nocc);

        //outfile->Printf(" Temp2:\n");
        //print_mat(Temp2,nb,nocc,outfile);

        double** Temp3 = block_matrix(na,nocc);
        C_DGEMM('N','N',na,nocc,nb,1.0,Sab[0],nb,Temp2[0],nocc,0.0,Temp3[0],nocc);

        //outfile->Printf(" Temp3:\n");
        //print_mat(Temp3,na,nocc,outfile);

        double** T = block_matrix(nocc,nocc);
        C_DGEMM('T','N',nocc,nocc,na,1.0,Ca[0],nocc,Temp3[0],nocc,0.0,T[0],nocc);

        //outfile->Printf(" T:\n");
        //print_mat(T,nocc,nocc,outfile);

        //Find T^-1/2
        // First, diagonalize T
        // the C_DSYEV call replaces the original matrix T with its eigenvectors
        double* eigval = init_array(nocc);
        int lwork = nocc * 3;
        double* work = init_array(lwork);
        int stat = C_DSYEV('v','u',nocc,T[0],nocc,eigval, work,lwork);
        if (stat != 0) {
            outfile->Printf( "C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }
        free(work);

        // Now T contains the eigenvectors of the original T
        // Copy T to T_copy
        double **T_mhalf = block_matrix(nocc, nocc);
        double **T_copy = block_matrix(nocc, nocc);
        C_DCOPY(nocc*nocc,T[0],1,T_copy[0],1);

        // Now form T^{-1/2} = U(T)*t^{-1/2}*U,
        // where t^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of T
        for (int i=0; i<nocc; i++) {
            if (eigval[i] < 1.0E-10)
                eigval[i] = 0.0;
            else
                eigval[i] = 1.0 / sqrt(eigval[i]);

            // scale one set of eigenvectors by the diagonal elements t^{-1/2}
            C_DSCAL(nocc, eigval[i], T[i], 1);
        }
        free(eigval);

        // T_mhalf = T_copy(T) * T
        C_DGEMM('t','n',nocc,nocc,nocc,1.0,
                T_copy[0],nocc,T[0],nocc,0.0,T_mhalf[0],nocc);

        //Form CB
        C_DGEMM('N','N',nb,nocc,nocc,1.0,Temp2[0],nocc,T_mhalf[0],nocc,0.0,Cb[0],nocc);

        free_block(Temp1);
        free_block(Temp2);
        free_block(Temp3);
        free_block(T);
        free_block(T_copy);
        free_block(T_mhalf);

    }
    return C_B;
}

}}
